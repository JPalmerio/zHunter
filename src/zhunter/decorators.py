import functools
import logging

log = logging.getLogger(__name__)


def check_active(func):
    """Check that the widget (passed in self) is active
    before running the function.
    """
    @functools.wraps(func)
    def wrapper_check_active(self, *args, **kwargs):
        if self.active:
            func(self, *args, **kwargs)
            return func(self, *args, **kwargs)
        else:
            log.debug(f"Ignoring call to {func.__name__} as {self} is not active.")

    return wrapper_check_active


def debug(func):
    """Print the function signature and return value"""
    @functools.wraps(func)
    def wrapper_debug(*args, **kwargs):
        args_repr = [repr(a) for a in args]                      # 1
        kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]  # 2
        signature = ", ".join(args_repr + kwargs_repr)           # 3
        log.debug(f"Calling {func.__name__}({signature})")
        value = func(*args, **kwargs)
        log.debug(f"{func.__name__!r} returned {value!r}")           # 4
        return value
    return wrapper_debug