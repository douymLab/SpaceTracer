
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
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
) -> List[Any]:
    """
    Args:
        items: 任务列表
        worker_fn: 单任务处理函数
        max_workers: 最大并发数
        desc: 日志描述
        raise_on_error: 是否在任务失败时抛异常
        backend: 并行后端，"thread" 或 "process"

    Returns:
        与 items 对应顺序一致的结果列表
    """
    items = list(items)
    if not items:
        return []

    if backend not in {"thread", "process"}:
        raise ValueError(f"Unsupported parallel backend: {backend}")

    results = [None] * len(items)
    future_to_index = {}

    executor_cls = ThreadPoolExecutor if backend == "thread" else ProcessPoolExecutor
    with executor_cls(max_workers=min(len(items), max_workers)) as executor:
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

    if errors and raise_on_error:
        msg = "; ".join([f"index={i}, item={item}, err={err}" for i, item, err in errors])
        raise RuntimeError(f"parallel_map failed: {msg}")

    return results
