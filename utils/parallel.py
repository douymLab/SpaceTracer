from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
import time
from typing import Callable, Iterable, List, Any, Optional, Literal, Tuple, Sequence, Union
import logging

from SpaceTracer.utils.logger import get_logger, ProgressLogger

Backend = Literal["thread", "process"]

def parallel_map(
    items: Iterable[Any],
    worker_fn: Callable[..., Any],
    max_workers: int = 4,
    desc: Optional[str] = None,
    raise_on_error: bool = True,
    backend: Backend = "thread",
    progress_interval: float = 0.05,
    logger: Optional[logging.Logger] = None,
    worker_takes_tuple: bool = False,
) -> List[Any]:
    """
    Args:
        items: 任务列表
            - 如果 worker_takes_tuple=False: item 会直接传给 worker_fn(item)
            - 如果 worker_takes_tuple=True: item 必须是 tuple/list，会调用 worker_fn(*item)
        worker_fn: 单任务处理函数（需要可 pickle 才能用于 process）
        max_workers: 最大并发数
        desc: 进度描述
        raise_on_error: 是否在任务失败时抛异常（True=最后统一抛；False=失败返回 None）
        backend: 并行后端，"thread" 或 "process"
        progress_interval: 进度打印间隔（0-1之间），0.1=每10%更新；None=不显示
        logger: SpaceTracer logger；不传则 get_logger("parallel")
        worker_takes_tuple: 是否把 item 作为 *args 展开传入 worker_fn
    """
    items = list(items)
    if not items:
        return []

    if backend not in {"thread", "process"}:
        raise ValueError(f"Unsupported parallel backend: {backend}")

    log = logger or get_logger("parallel")
    label = desc or "parallel"

    total = len(items)
    results: List[Any] = [None] * total
    future_to_index = {}

    completed = 0
    start_time = time.time()

    use_progress = (progress_interval is not None)
    interval_pct = int(progress_interval * 100) if use_progress else None

    if use_progress:
        try:
            log.info(f"[{label}] Starting {total} tasks... (backend={backend}, max_workers={min(total, max_workers)})")
        except Exception:
            pass

    # 可选：如果你仍想保留 ProgressLogger 实例用于别处，这里留着
    prog = ProgressLogger(total=total, desc=f"[{label}]") if use_progress else None
    _ = prog  # 避免 lint 报未使用；我们用自定义 stderr bar 控制粒度

    def _submit(executor, item):
        if worker_takes_tuple:
            if not isinstance(item, (tuple, list)):
                raise TypeError(f"worker_takes_tuple=True but item is not tuple/list: {type(item)}")
            return executor.submit(worker_fn, *item)
        else:
            return executor.submit(worker_fn, item)

    executor_cls = ThreadPoolExecutor if backend == "thread" else ProcessPoolExecutor

    errors: List[Tuple[int, Any, Exception]] = []
    last_print_pct = -1

    with executor_cls(max_workers=min(total, max_workers)) as executor:
        for i, item in enumerate(items):
            future = _submit(executor, item)
            future_to_index[future] = i

        for future in as_completed(future_to_index):
            i = future_to_index[future]
            try:
                results[i] = future.result()
            except Exception as e:
                log.exception(f"[{label}] task failed at index={i}, item={items[i]}")
                errors.append((i, items[i], e))
                results[i] = None  # 失败时占位

            completed += 1

            if use_progress:
                current_pct = int(completed / total * 100)
                if interval_pct is None:
                    should_print = True
                else:
                    should_print = (current_pct == 100) or (current_pct >= last_print_pct + interval_pct)

                if should_print:
                    elapsed = time.time() - start_time
                    eta = (elapsed / completed) * (total - completed) if completed > 0 else 0.0

                    import sys
                    bar_len = 20
                    filled = int(bar_len * completed / total)
                    bar = '█' * filled + '░' * (bar_len - filled)
                    sys.stderr.write(f"[{label}] {bar} {current_pct}% ({completed}/{total}) | ETA: {eta:.0f}s\n")
                    sys.stderr.flush()

                    last_print_pct = current_pct

    if use_progress:
        import sys
        sys.stderr.write("\n")
        sys.stderr.flush()
    
    if errors and raise_on_error:
        msg = "; ".join([f"index={i},  err={err}" for i, item, err in errors])
        raise RuntimeError(f"parallel_map failed: {msg}")

    if use_progress:
        elapsed = time.time() - start_time
        log.info(f"[{label}] Completed! {total} tasks in {elapsed:.1f}s")

    return results
