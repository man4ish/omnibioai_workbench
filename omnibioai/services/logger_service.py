import logging

logger = logging.getLogger("OmniBioAI")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# Optional utility functions
def ensure_dir(path):
    import os
    if not os.path.exists(path):
        os.makedirs(path)

