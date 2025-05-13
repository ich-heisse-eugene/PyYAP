import logging
import logging.handlers

def listener(queue, logfile, level=logging.INFO):
    logger = logging.getLogger()
    logger.setLevel(level)
    fh = logging.FileHandler(logfile, "a")
    fmt = logging.Formatter('%(asctime)s - %(message)s')
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    while True:
        record = queue.get()
        if record is None:
            break
        logger.log(record[0], record[1])
    return None
