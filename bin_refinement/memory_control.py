import tracemalloc
import logging

def sizeof_fmt(num, suffix="B"):
    for unit in ["", "K", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"

def measure_memory(function):

    def new_function(*args, **kwargs):

        tracemalloc.start()
        ret = function(*args, **kwargs)
        current, peak = tracemalloc.get_traced_memory()

        logging.info(f'In  function {function.__name__}')
        logging.info(f"\033[37mCurrent memory usage:\033[36m {sizeof_fmt(current)}\033[0m")
        logging.info(f"\033[37mPeak                :\033[36m {sizeof_fmt(peak)}\033[0m")
        tracemalloc.stop()
        return ret

    return new_function



