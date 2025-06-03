# src/palantir/logger.py

import logging

def setup_logging(log_filepath: str = "", verbose: bool = False):
    logger = logging.getLogger("palantir")
    logger.setLevel(logging.DEBUG)

    if logger.handlers:
        return 

    file_handler = logging.FileHandler(log_filepath + "palantir_info.log")
    file_handler.setLevel(logging.DEBUG)
    file_format = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    file_handler.setFormatter(file_format)
    logger.addHandler(file_handler)

    if verbose:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_format = logging.Formatter('%(levelname)s - %(message)s')
        console_handler.setFormatter(console_format)
        logger.addHandler(console_handler)
