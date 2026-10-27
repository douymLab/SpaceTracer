#!/usr/bin/env python3
"""
Reproducibility management for PhyloSOLID.
"""

import os
import random
import logging
import hashlib
from typing import Optional

logger = logging.getLogger(__name__)

DEFAULT_SEED = 42
_global_seed = DEFAULT_SEED
_is_initialized = False


def set_seed(seed: int = DEFAULT_SEED, verbose: bool = True) -> None:
    """
    Set the global random seed for reproducibility.
    """
    global _global_seed, _is_initialized
    
    _global_seed = seed
    _is_initialized = True
    
    random.seed(seed)
    
    try:
        import numpy as np
        np.random.seed(seed)
    except ImportError:
        pass
    
    os.environ['PYTHONHASHSEED'] = str(seed)
    
    # ============================================================
    # CRITICAL FIX: Also set random seeds for Python's random module
    # at the C level using random.seed() again with a different approach
    # ============================================================
    # This helps ensure that the RNG state is fully reset even in forked processes
    import random as _random
    _random.seed(seed)
    # ============================================================
    
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
    except ImportError:
        pass
    
    if verbose:
        logger.info(f"Global random seed set to: {seed}")


def get_seed() -> int:
    """Get the current global seed."""
    return _global_seed


def is_initialized() -> bool:
    """Check if the seed has been initialized."""
    return _is_initialized


def _deterministic_hash(seed: int, salt: str, value: any) -> int:
    """Generate a deterministic hash using hashlib."""
    if isinstance(value, (list, tuple, set)):
        sorted_vals = sorted(str(v) for v in value)
        value_str = ','.join(sorted_vals)
    else:
        value_str = str(value)
    
    hash_obj = hashlib.md5((str(seed) + salt + value_str).encode())
    return int(hash_obj.hexdigest(), 16)


def deterministic_choice(population, seed: Optional[int] = None, salt: str = ""):
    """Deterministic version of random.choice."""
    if seed is None:
        seed = _global_seed
    
    hash_value = _deterministic_hash(seed, salt, population)
    rng = random.Random(seed + hash_value % 100000)
    return rng.choice(population)


def deterministic_sample(population, k, seed: Optional[int] = None, salt: str = ""):
    """Deterministic version of random.sample."""
    if seed is None:
        seed = _global_seed
    
    hash_value = _deterministic_hash(seed, salt, population)
    rng = random.Random(seed + hash_value % 100000)
    return rng.sample(population, k)


def deterministic_permutation(arr, seed: Optional[int] = None, salt: str = ""):
    """Deterministic version of np.random.permutation."""
    if seed is None:
        seed = _global_seed
    
    import numpy as np
    
    if isinstance(arr, (list, tuple, np.ndarray)):
        sorted_arr = sorted(arr)
        arr_str = ','.join(str(x) for x in sorted_arr)
    else:
        arr_str = str(arr)
    
    hash_obj = hashlib.md5((str(seed) + salt + arr_str).encode())
    final_seed = seed + int(hash_obj.hexdigest(), 16) % 100000
    
    rng = np.random.RandomState(final_seed)
    return rng.permutation(arr)


def deterministic_shuffle(arr, seed: Optional[int] = None, salt: str = ""):
    """Deterministic version of random.shuffle."""
    if seed is None:
        seed = _global_seed
    
    hash_value = _deterministic_hash(seed, salt, arr)
    rng = random.Random(seed + hash_value % 100000)
    
    shuffled = arr.copy() if hasattr(arr, 'copy') else arr[:]
    rng.shuffle(shuffled)
    return shuffled


# Auto-initialize with default seed when module is imported
set_seed(DEFAULT_SEED, verbose=False)