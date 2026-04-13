
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
from typing import Callable, Iterable, List, Any, Optional
import logging

logger = logging.getLogger(__name__)

def parallel_map(
    items: Iterable[Any],
    worker_fn: Callable[[Any], Any],
    max_workers: int = 4,
    desc: Optional[str] = None,
    raise_on_error: bool = True,
    backend: str = "thread",
    progress_interval: float = 0.05,  # 新增：每 10% 打印一次
) -> List[Any]:
    """
    Args:
        items: 任务列表
        worker_fn: 单任务处理函数
        max_workers: 最大并发数
        desc: 日志描述
        raise_on_error: 是否在任务失败时抛异常
        backend: 并行后端，"thread" 或 "process"
        progress_interval: 进度打印间隔（0-1之间），0.1=每10%打印，None=不打印
    """
    items = list(items)
    if not items:
        return []

    if backend not in {"thread", "process"}:
        raise ValueError(f"Unsupported parallel backend: {backend}")

    total = len(items)
    results = [None] * total
    future_to_index = {}
    
    # 进度追踪
    completed = 0
    last_print_pct = 0
    start_time = time.time()
    
    if desc and progress_interval is not None:
        print(f"[{desc}] Starting {total} tasks...")

    executor_cls = ThreadPoolExecutor if backend == "thread" else ProcessPoolExecutor
    with executor_cls(max_workers=min(total, max_workers)) as executor:
        for i, item in enumerate(items):
            future = executor.submit(worker_fn, item)
            future_to_index[future] = i

        errors = []
        for future in as_completed(future_to_index):
            i = future_to_index[future]
            try:
                results[i] = future.result()
            except Exception as e:
                logger.exception(f"[parallel] task failed at index={i}, item={items[i]}")
                errors.append((i, items[i], e))
            
            # 进度更新
            completed += 1
            if progress_interval is not None:
                current_pct = int(completed / total * 100)
                interval_pct = int(progress_interval * 100)
                
                if current_pct >= last_print_pct + interval_pct or completed == total:
                    elapsed = time.time() - start_time
                    eta = (elapsed / completed) * (total - completed) if completed > 0 else 0
                    bar_len = 20
                    filled = int(bar_len * completed / total)
                    bar = '█' * filled + '░' * (bar_len - filled)
                    print(f"  [{desc}] {bar} {current_pct}% ({completed}/{total}) | Estimated Time Remaining: {eta:.0f}s")
                    last_print_pct = current_pct    


    if errors and raise_on_error:
        msg = "; ".join([f"index={i}, item={item}, err={err}" for i, item, err in errors])
        raise RuntimeError(f"parallel_map failed: {msg}")
    
    if desc and progress_interval is not None:
        elapsed = time.time() - start_time
        print(f"[{desc}] Completed! {total} tasks in {elapsed:.1f}s")

    return results