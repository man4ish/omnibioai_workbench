"""
base.py

Abstract base class for workflow engine adapters.
"""

from abc import ABC, abstractmethod
from typing import Iterator, Tuple


class BaseWorkflowAdapter(ABC):
    """
    Base adapter contract for workflow engines.
    """

    @abstractmethod
    def run(
        self,
        entrypoint: str,
        work_dir: str,
        params: dict
    ) -> Iterator[Tuple[int, str]]:
        """
        Execute the workflow.

        Yields:
            (progress_percent, log_line)

        Notes:
            Progress may be heuristic.
        """
        pass
