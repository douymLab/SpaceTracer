from concurrent.futures import (
    ThreadPoolExecutor,
    ProcessPoolExecutor,
    wait,
    FIRST_COMPLETED,
)
import time
from typing import Callable, Iterable, List, Any, Optional, Literal, Tuple, Generator
import logging
import os
import sys
import gc

import psutil

from SpaceTracer.utils.logger import get_logger, ProgressLogger

Backend = Literal["thread", "process"]


def _safe_mem_info_mb(proc: psutil.Process) -> dict:
    """
    返回进程内存信息（MB）:
    - rss: 常驻内存
    - uss: 独占内存（如果系统支持）
    """
    rss_mb = 0.0
    uss_mb = None

    try:
        rss_mb = proc.memory_info().rss / 1024 / 1024
    except Exception:
        pass

    try:
        full = proc.memory_full_info()
        uss = getattr(full, "uss", None)
        if uss is not None:
            uss_mb = uss / 1024 / 1024
    except Exception:
        pass

    return {
        "rss_mb": rss_mb,
        "uss_mb": uss_mb,
    }


def _get_process_mem_summary(pid: Optional[int] = None) -> dict:
    """
    获取单个进程内存摘要
    """
    pid = pid or os.getpid()
    proc = psutil.Process(pid)
    mem = _safe_mem_info_mb(proc)
    return {
        "pid": pid,
        "rss_mb": mem["rss_mb"],
        "uss_mb": mem["uss_mb"],
    }


def _get_children_mem_summary(top_n: int = 5) -> dict:
    """
    获取当前主进程所有子进程的内存摘要
    """
    parent = psutil.Process(os.getpid())

    total_rss_mb = 0.0
    total_uss_mb = 0.0
    uss_available = True

    children_rows = []

    for child in parent.children(recursive=True):
        try:
            mem = _safe_mem_info_mb(child)
            rss_mb = mem["rss_mb"]
            uss_mb = mem["uss_mb"]

            total_rss_mb += rss_mb
            if uss_mb is None:
                uss_available = False
            else:
                total_uss_mb += uss_mb

            try:
                name = child.name()
            except Exception:
                name = "unknown"

            try:
                cmdline = " ".join(child.cmdline()[:8])
            except Exception:
                cmdline = ""

            children_rows.append({
                "pid": child.pid,
                "name": name,
                "rss_mb": rss_mb,
                "uss_mb": uss_mb,
                "cmdline": cmdline,
            })
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue

    children_rows.sort(key=lambda x: x["rss_mb"], reverse=True)

    return {
        "n_children": len(children_rows),
        "children_rss_mb": total_rss_mb,
        "children_uss_mb": (total_uss_mb if uss_available else None),
        "top_children": children_rows[:top_n],
    }


def _format_top_children(top_children: List[dict], use_uss: bool = False) -> str:
    """
    把 top children 格式化成一行短文本
    """
    if not top_children:
        return "None"

    parts = []
    for row in top_children:
        mem_val = row["uss_mb"] if use_uss and row["uss_mb"] is not None else row["rss_mb"]
        mem_key = "USS" if use_uss and row["uss_mb"] is not None else "RSS"
        parts.append(f"{row['name']}:{row['pid']}({mem_key}={mem_val:.1f}MB)")
    return ", ".join(parts)


def get_memory_snapshot(top_n: int = 5) -> dict:
    """
    获取当前主进程 + 子进程的完整内存快照
    """
    main = _get_process_mem_summary()
    children = _get_children_mem_summary(top_n=top_n)

    total_rss_mb = main["rss_mb"] + children["children_rss_mb"]

    total_uss_mb = None
    if main["uss_mb"] is not None and children["children_uss_mb"] is not None:
        total_uss_mb = main["uss_mb"] + children["children_uss_mb"]

    return {
        "main_pid": main["pid"],
        "main_rss_mb": main["rss_mb"],
        "main_uss_mb": main["uss_mb"],
        "children_rss_mb": children["children_rss_mb"],
        "children_uss_mb": children["children_uss_mb"],
        "total_rss_mb": total_rss_mb,
        "total_uss_mb": total_uss_mb,
        "n_children": children["n_children"],
        "top_children": children["top_children"],
    }


def memory_checkpoint(
    label: str,
    logger: Optional[logging.Logger] = None,
    top_n: int = 5,
    print_children: bool = True,
    do_gc: bool = False,
) -> None:
    """
    在任意 step 前后打一个内存检查点
    """
    log = logger or get_logger("memory")

    if do_gc:
        gc.collect()

    snap = get_memory_snapshot(top_n=top_n)

    # main_uss_text = f"{snap['main_uss_mb']:.1f}MB" if snap["main_uss_mb"] is not None else "N/A"
    # child_uss_text = f"{snap['children_uss_mb']:.1f}MB" if snap["children_uss_mb"] is not None else "N/A"
    total_uss_text = f"{snap['total_uss_mb']:.1f}MB" if snap["total_uss_mb"] is not None else "N/A"

    log.info(
        f"[MEM][{label}] "
        # f"Main(pid={snap['main_pid']}) RSS={snap['main_rss_mb']:.1f}MB USS={main_uss_text} | "
        # f"Children RSS={snap['children_rss_mb']:.1f}MB USS={child_uss_text} n={snap['n_children']} | "
        f"Total RSS={snap['total_rss_mb']:.1f}MB USS={total_uss_text}"
    )

    # if print_children and snap["top_children"]:
    #     top_rss = _format_top_children(snap["top_children"], use_uss=False)
    #     log.info(f"[MEM][{label}] Top children by RSS: {top_rss}")

    #     has_any_uss = any(row["uss_mb"] is not None for row in snap["top_children"])
    #     if has_any_uss:
    #         top_uss = _format_top_children(snap["top_children"], use_uss=True)
    #         log.info(f"[MEM][{label}] Top children by USS: {top_uss}")


def _read_int_file(path: str) -> Optional[int]:
    try:
        with open(path, "r") as f:
            raw = f.read().strip()
        if raw in {"", "max"}:
            return None
        return int(raw)
    except Exception:
        return None


def _get_memory_usage_bytes() -> Optional[int]:
    """
    获取当前总内存使用量（bytes）
    优先 cgroup，其次 fallback 到 main + children RSS
    """
    usage = _read_int_file("/sys/fs/cgroup/memory.current")
    if usage is not None:
        return usage

    usage = _read_int_file("/sys/fs/cgroup/memory/memory.usage_in_bytes")
    if usage is not None:
        return usage

    try:
        snap = get_memory_snapshot(top_n=0)
        return int(snap["total_rss_mb"] * 1024 * 1024)
    except Exception:
        return None


def _get_memory_limit_bytes() -> Optional[int]:
    """
    获取 cgroup memory limit（bytes）
    """
    limit = _read_int_file("/sys/fs/cgroup/memory.max")
    if limit is not None and limit > 0:
        return limit

    limit = _read_int_file("/sys/fs/cgroup/memory/memory.limit_in_bytes")
    if limit is not None and 0 < limit < (1 << 60):
        return limit

    return None


def _worker_wrapper(worker_fn, item, worker_takes_tuple: bool):
    """
    worker 包装器：
    - 执行任务
    - 返回任务结果
    - 顺便记录 worker 进程内存信息
    """
    pid = os.getpid()

    if worker_takes_tuple:
        if not isinstance(item, (tuple, list)):
            raise TypeError(f"worker_takes_tuple=True but item is not tuple/list: {type(item)}")
        result = worker_fn(*item)
    else:
        result = worker_fn(item)

    mem = _get_process_mem_summary(pid)

    return {
        "result": result,
        "worker_pid": pid,
        "worker_rss_mb": mem["rss_mb"],
        "worker_uss_mb": mem["uss_mb"],
    }


def _parallel_core(
    items: Iterable[Any],
    worker_fn: Callable[..., Any],
    max_workers: int = 4,
    desc: Optional[str] = None,
    raise_on_error: bool = True,
    backend: Backend = "thread",
    progress_interval: float = 0.05,
    logger: Optional[logging.Logger] = None,
    worker_takes_tuple: bool = False,
    top_n_children: int = 3,
    debug_children_every_print: bool = False,
    max_in_flight: Optional[int] = None,
    memory_limit_bytes=None,
    memory_soft_ratio=0.85,
    memory_check_interval=0.5,
    wait_on_memory_pressure=True,
) -> Generator[Tuple[int, Any], None, None]:
    """
    内部统一调度核心
    流式 yield: (index, result)
    """
    items = list(items)
    if not items:
        return

    if backend not in {"thread", "process"}:
        raise ValueError(f"Unsupported parallel backend: {backend}")

    log = logger or get_logger("parallel")
    label = desc or "parallel"

    total = len(items)
    completed = 0
    start_time = time.time()

    use_progress = (progress_interval is not None)
    interval_pct = int(progress_interval * 100) if use_progress else None

    if use_progress:
        try:
            log.info(
                f"[{label}] Starting {total} tasks... "
                f"(backend={backend}, max_workers={min(total, max_workers)})"
            )
        except Exception:
            pass

    prog = ProgressLogger(total=total, desc=f"[{label}]") if use_progress else None
    _ = prog

    executor_cls = ThreadPoolExecutor if backend == "thread" else ProcessPoolExecutor

    errors: List[Tuple[int, Any, Exception]] = []
    last_print_pct = -1

    n_workers = min(total, max_workers)
    if max_in_flight is None:
        max_in_flight = total
    else:
        max_in_flight = max(n_workers, min(max_in_flight, total))

    if wait_on_memory_pressure and memory_limit_bytes is None:
        memory_limit_bytes = _get_memory_limit_bytes()

    soft_limit_bytes = None
    if wait_on_memory_pressure and memory_limit_bytes is not None:
        soft_limit_bytes = int(memory_limit_bytes * memory_soft_ratio)
        log.info(
            f"[{label}] memory-aware submit enabled: "
            f"limit={memory_limit_bytes / 1024**3:.2f}GB, "
            f"soft={soft_limit_bytes / 1024**3:.2f}GB"
        )

    last_memory_check_time = 0.0
    last_memory_usage_bytes = None
    last_pressure_log_time = 0.0

    def _memory_allows_submit() -> bool:
        nonlocal last_memory_check_time, last_memory_usage_bytes

        if not wait_on_memory_pressure or soft_limit_bytes is None:
            return True

        now = time.time()
        if (
            last_memory_usage_bytes is None
            or now - last_memory_check_time >= memory_check_interval
        ):
            # 修正：这里必须是 bytes
            last_memory_usage_bytes = _get_memory_usage_bytes()
            last_memory_check_time = now

        if last_memory_usage_bytes is None:
            return True

        return last_memory_usage_bytes < soft_limit_bytes

    def _print_progress(task_idx, task_status, worker_pid, worker_rss_mb, worker_uss_mb):
        nonlocal last_print_pct

        if not use_progress:
            return

        current_pct = int(completed / total * 100)
        if interval_pct is None:
            should_print = True
        else:
            should_print = (current_pct == 100) or (current_pct >= last_print_pct + interval_pct)

        if not should_print:
            return

        elapsed = time.time() - start_time
        eta = (elapsed / completed) * (total - completed) if completed > 0 else 0.0
        current_time = time.strftime("%H:%M:%S", time.localtime())

        bar_len = 20
        filled = int(bar_len * completed / total)
        bar = '█' * filled + '░' * (bar_len - filled)

        snap = get_memory_snapshot(top_n=top_n_children)

        main_uss_text = f"{snap['main_uss_mb']:.1f}MB" if snap["main_uss_mb"] is not None else "N/A"
        child_uss_text = f"{snap['children_uss_mb']:.1f}MB" if snap["children_uss_mb"] is not None else "N/A"
        total_uss_text = f"{snap['total_uss_mb']:.1f}MB" if snap["total_uss_mb"] is not None else "N/A"

        worker_rss_text = f"{worker_rss_mb:.1f}MB" if worker_rss_mb is not None else "N/A"
        worker_uss_text = f"{worker_uss_mb:.1f}MB" if worker_uss_mb is not None else "N/A"

        # sys.stderr.write(
        #     f"[{label}] {bar} {current_pct}% ({completed}/{total}) | "
        #     f"task_id={task_idx} | status={task_status} | worker_pid={worker_pid} | "
        #     f"workerRSS={worker_rss_text} workerUSS={worker_uss_text} | "
        #     f"MainRSS={snap['main_rss_mb']:.1f}MB MainUSS={main_uss_text} | "
        #     f"ChildrenRSS={snap['children_rss_mb']:.1f}MB ChildrenUSS={child_uss_text} | "
        #     f"TotalRSS={snap['total_rss_mb']:.1f}MB TotalUSS={total_uss_text} | "
        #     f"n_children={snap['n_children']} | "
        #     f"in_flight<={max_in_flight} | "
        #     f"Cost: {elapsed:.1f}s | ETA: {eta:.1f}s | {current_time}\n"
        # )
        sys.stderr.write(
            f"[{label}] {bar} {current_pct}% ({completed}/{total}) | "
            f"task_id={task_idx} | status={task_status} | "
            # f"workerRSS={worker_rss_text} workerUSS={worker_uss_text} | "
            # f"MainRSS={snap['main_rss_mb']:.1f}MB MainUSS={main_uss_text} | "
            # f"ChildrenRSS={snap['children_rss_mb']:.1f}MB ChildrenUSS={child_uss_text} | "
            f"TotalRSS={snap['total_rss_mb']:.1f}MB TotalUSS={total_uss_text} | "
            # f"n_children={snap['n_children']} | "
            f"runing_tasks<={max_in_flight} | "
            f"Cost: {elapsed:.1f}s | ETA: {eta:.1f}s | {current_time}\n"
        )
        sys.stderr.flush()

        if debug_children_every_print and snap["top_children"]:
            top_rss = _format_top_children(snap["top_children"], use_uss=False)
            sys.stderr.write(f"[{label}] top_children_rss: {top_rss}\n")
            sys.stderr.flush()

        last_print_pct = current_pct

    with executor_cls(max_workers=n_workers) as executor:
        future_to_index = {}
        next_submit_idx = 0

        def submit_one(idx: int):
            item = items[idx]
            future = executor.submit(_worker_wrapper, worker_fn, item, worker_takes_tuple)
            future_to_index[future] = idx

        def try_submit_one():
            nonlocal next_submit_idx, last_pressure_log_time
            if next_submit_idx >= total:
                return False
            if len(future_to_index) >= max_in_flight:
                return False

            if future_to_index and not _memory_allows_submit():
                now = time.time()
                if now - last_pressure_log_time >= 2.0:
                    if last_memory_usage_bytes is not None and soft_limit_bytes is not None:
                        log.info(
                            f"[{label}] memory pressure high, pause submit: "
                            f"used={last_memory_usage_bytes / 1024**3:.2f}GB, "
                            f"soft={soft_limit_bytes / 1024**3:.2f}GB, "
                            f"in_flight={len(future_to_index)}"
                        )
                    last_pressure_log_time = now
                return False

            submit_one(next_submit_idx)
            next_submit_idx += 1
            return True

        while try_submit_one():
            pass

        while future_to_index:
            done, _ = wait(future_to_index.keys(), return_when=FIRST_COMPLETED)

            for future in done:
                i = future_to_index.pop(future)

                worker_pid = None
                worker_rss_mb = None
                worker_uss_mb = None
                task_status = "OK"
                result = None

                try:
                    payload = future.result()
                    result = payload["result"]
                    worker_pid = payload.get("worker_pid")
                    worker_rss_mb = payload.get("worker_rss_mb")
                    worker_uss_mb = payload.get("worker_uss_mb")
                except Exception as e:
                    log.exception(f"[{label}] task failed at index={i}, item={items[i]}")
                    errors.append((i, items[i], e))
                    task_status = "FAIL"

                completed += 1
                _print_progress(i, task_status, worker_pid, worker_rss_mb, worker_uss_mb)

                yield i, result

            last_memory_usage_bytes = _get_memory_usage_bytes()
            last_memory_check_time = time.time()

            while try_submit_one():
                pass

    if use_progress:
        sys.stderr.write("\n")
        sys.stderr.flush()

    if errors and raise_on_error:
        msg = "; ".join([f"index={i}, err={err}" for i, item, err in errors])
        raise RuntimeError(f"{label} failed: {msg}")

    if use_progress:
        elapsed = time.time() - start_time
        log.info(f"[{label}] Completed! {completed} tasks in {elapsed:.1f}s")


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
    top_n_children: int = 3,
    debug_children_every_print: bool = False,
    max_in_flight: Optional[int] = None,
    memory_limit_bytes=None,
    memory_soft_ratio=0.85,
    memory_check_interval=0.5,
    wait_on_memory_pressure=True,
) -> List[Any]:
    items = list(items)
    if not items:
        return []

    results: List[Any] = [None] * len(items)

    for i, result in _parallel_core(
        items=items,
        worker_fn=worker_fn,
        max_workers=max_workers,
        desc=desc,
        raise_on_error=raise_on_error,
        backend=backend,
        progress_interval=progress_interval,
        logger=logger,
        worker_takes_tuple=worker_takes_tuple,
        top_n_children=top_n_children,
        debug_children_every_print=debug_children_every_print,
        max_in_flight=max_in_flight,
        memory_limit_bytes=memory_limit_bytes,
        memory_soft_ratio=memory_soft_ratio,
        memory_check_interval=memory_check_interval,
        wait_on_memory_pressure=wait_on_memory_pressure,
    ):
        results[i] = result

    return results


def parallel_imap(
    items: Iterable[Any],
    worker_fn: Callable[..., Any],
    max_workers: int = 4,
    desc: Optional[str] = None,
    raise_on_error: bool = True,
    backend: Backend = "thread",
    progress_interval: float = 0.05,
    logger: Optional[logging.Logger] = None,
    worker_takes_tuple: bool = False,
    top_n_children: int = 3,
    debug_children_every_print: bool = False,
    max_in_flight: Optional[int] = None,
    memory_limit_bytes=None,
    memory_soft_ratio=0.85,
    memory_check_interval=0.5,
    wait_on_memory_pressure=True,
) -> Generator[Tuple[int, Any], None, None]:

    yield from _parallel_core(
        items=items,
        worker_fn=worker_fn,
        max_workers=max_workers,
        desc=desc,
        raise_on_error=raise_on_error,
        backend=backend,
        progress_interval=progress_interval,
        logger=logger,
        worker_takes_tuple=worker_takes_tuple,
        top_n_children=top_n_children,
        debug_children_every_print=debug_children_every_print,
        max_in_flight=max_in_flight,
        memory_limit_bytes=memory_limit_bytes,
        memory_soft_ratio=memory_soft_ratio,
        memory_check_interval=memory_check_interval,
        wait_on_memory_pressure=wait_on_memory_pressure,
    )
