
# SpaceTracer/utils/logger.py
"""
Logger setup for SpaceTracer
"""

import logging
import sys
from typing import Optional
from pathlib import Path
from colorama import init, Fore, Style

# 定义颜色（可选）
# try:
#     init(autoreset=True)
#     COLORS = {
#         'DEBUG': Fore.CYAN,
#         'INFO': Fore.GREEN,
#         'WARNING': Fore.YELLOW,
#         'ERROR': Fore.RED,
#         'CRITICAL': Fore.RED + Style.BRIGHT,
#     }
# except ImportError:
COLORS = {}


class ColoredFormatter(logging.Formatter):
    """logging formatter with color"""
    
    def __init__(self, fmt=None, datefmt=None, style='%'):
        super().__init__(fmt, datefmt, style)
    
    def format(self, record):
        """add color"""
        original_levelname = record.levelname # keep original log level
        if COLORS and record.levelname in COLORS: # add color
            record.levelname = f"{COLORS[record.levelname]}{record.levelname}{Style.RESET_ALL}"
            record.msg = f"{COLORS[record.levelname]}{record.msg}{Style.RESET_ALL}"
        
        result = super().format(record) # use the parent format of log
        record.levelname = original_levelname # but keep the original levelname of log
        
        return result


class ProgressFilter(logging.Filter):
    """handle progress item"""
    
    def filter(self, record):
        # if the log info stared with \n, which represents the progress item, the log info would not be record
        if record.getMessage().startswith('\r'):
            return False
        return True


def setup_logger(
    level: str = 'INFO',
    log_file: Optional[str] = None,
    verbose: bool = False,
    debug: bool = False
) -> logging.Logger:
    """
    set the logger for the program
    
    Args:
        level: log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: the log file path. If None, the log will be printed
        verbose: do you want to export detailed log info? 
        debug: do you want to use debug mode?     
    Returns:
        logger
    """
    # 确定日志级别
    if debug:
        log_level = logging.DEBUG
    elif verbose:
        log_level = logging.INFO
    else:
        log_level = getattr(logging, level.upper(), logging.WARNING)
    
    # 获取根logger
    logger = logging.getLogger('SpaceTracer')
    logger.setLevel(log_level)
    
    # 避免重复添加handler
    if logger.handlers:
        return logger
    
    # 控制台handler
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(log_level)
    
    # 根据级别设置不同的控制台格式
    if log_level <= logging.DEBUG:
        # 调试模式：显示详细信息
        console_format = ColoredFormatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s',
            datefmt='%H:%M:%S'
        )
    elif log_level <= logging.INFO:
        # 信息模式：显示基本信息
        console_format = ColoredFormatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%H:%M:%S'
        )
    else:
        # 警告及以上：只显示消息
        console_format = ColoredFormatter(
            '%(levelname)s: %(message)s'
        )
    
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # 文件handler（如果指定了日志文件）
    if log_file:
        # 确保日志目录存在
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)  # 文件记录所有级别的日志
        
        # 文件格式总是完整的
        file_format = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_format)
        file_handler.addFilter(ProgressFilter())  # 过滤进度条
        logger.addHandler(file_handler)
    
    # 避免传播到根logger（防止重复输出）
    logger.propagate = False
    
    # 记录初始化信息
    if debug:
        logger.debug(f"Logger initialized with level: {level}")
        logger.debug(f"Log file: {log_file if log_file else 'Not set'}")
    
    return logger


def get_logger(name: str = None) -> logging.Logger:
    """
    获取指定名称的logger
    
    Args:
        name: logger名称，如果为None则返回根logger
    
    Returns:
        logger实例
    """
    if name:
        return logging.getLogger(f'SpaceTracer.{name}')
    return logging.getLogger('SpaceTracer')


# 进度条相关的日志工具
class ProgressLogger:
    """进度条记录器"""
    
    def __init__(self, total: int, desc: str = "Processing", logger=None):
        self.total = total
        self.desc = desc
        self.current = 0
        self.logger = logger or get_logger('progress')
        self._last_percentage = -1
    
    def update(self, increment: int = 1):
        """更新进度"""
        self.current += increment
        percentage = int((self.current / self.total) * 100)

        # 每5%或完成时更新
        if percentage != self._last_percentage and (percentage % 5 == 0 or self.current == self.total):
            self._last_percentage = percentage
            if self.current == self.total:
                self.logger.info(f"{self.desc}: 100% ({self.current}/{self.total}) - Completed")
            else:
                # 使用\r实现进度条效果
                sys.stderr.write(f"{self.desc}: {percentage}% ({self.current}/{self.total})")
                sys.stderr.flush()
    
    def finish(self):
        """完成进度条"""
        if self.current < self.total:
            self.current = self.total
            self.logger.info(f"{self.desc}: 100% ({self.current}/{self.total}) - Completed")
        sys.stderr.write('\n')
        sys.stderr.flush()


# 装饰器：记录函数执行时间
def log_execution_time(logger=None):
    """
    记录函数执行时间的装饰器
    
    Example:
        @log_execution_time()
        def my_function():
            pass
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            import time
            start_time = time.time()
            
            log = logger or get_logger(func.__module__)
            log.debug(f"Starting {func.__name__}")
            
            try:
                result = func(*args, **kwargs)
                elapsed = time.time() - start_time
                log.debug(f"Finished {func.__name__} in {elapsed:.2f}s")
                return result
            except Exception as e:
                elapsed = time.time() - start_time
                log.error(f"{func.__name__} failed after {elapsed:.2f}s: {e}")
                raise
        
        return wrapper
    return decorator

