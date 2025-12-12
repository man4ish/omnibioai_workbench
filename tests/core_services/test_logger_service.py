from omnibioai.services.logger_service import logger, ensure_dir
import os

# ----------------------------
# 1. Test logging output
# ----------------------------
logger.info("This is an INFO test message")
logger.warning("This is a WARNING test message")
logger.error("This is an ERROR test message")

# ----------------------------
# 2. Test ensure_dir utility
# ----------------------------
test_dir = "data/test_logger_dir"
ensure_dir(test_dir)

if os.path.exists(test_dir) and os.path.isdir(test_dir):
    logger.info(f"ensure_dir successfully created directory: {test_dir}")
else:
    logger.error(f"ensure_dir failed to create directory: {test_dir}")

# ----------------------------
# 3. Cleanup (optional)
# ----------------------------
import shutil
shutil.rmtree(test_dir)
logger.info(f"Cleaned up test directory: {test_dir}")

