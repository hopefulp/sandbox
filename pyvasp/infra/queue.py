from dataclasses import dataclass
from typing import Optional

@dataclass
class QueueConfig:
    partition: int
    nnode: int
    nproc: Optional[int] = None
