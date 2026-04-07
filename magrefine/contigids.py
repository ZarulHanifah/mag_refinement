
from __future__ import annotations

import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class BaseContigID(ABC):
    """
    Abstract base class representing a single contig, parsing its metadata
    from the FASTA header and holding its abundance information.
    """
    _name: str  
    _abund_info: dict[str, float] = field(default_factory=dict)

    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    @abstractmethod
    def length(self) -> int | None: ...

    def get_abund_info(self) -> dict[str, float]:
        """Returns the abundance information for this contig."""
        return self._abund_info

    @property
    def depth_from_all_samples(self) -> float:
        if not self._abund_info:
            return 0.0
        return sum(self._abund_info.values())

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name={self.name}, length={self.length})"



@dataclass(frozen=True)
class MyloParsedHeader:
    """A struct to hold the parsed data from a myloasm contig header."""
    name: str
    length: int
    is_circular: bool
    depth: str
    is_duplicated: bool


@dataclass
class MyloContigID(BaseContigID):
    """
    Represents a single contig from myloasm, parsing its specific rich metadata 
    from the FASTA header.
    """
    length_val: int = 0
    is_circular: bool = False
    mylo_depth: str = ""
    is_duplicated: bool = False

    HEADER_PATTERN = re.compile(
        r"^>?(?P<name>\w+)_len-(?P<length>\d+)_circular-(?P<circular>yes|possibly|no)_depth-(?P<depth>[\d.-]+)_duplicated-(?P<duplicated>yes|possibly|no)$"
    )

    @classmethod
    def from_header(cls, header: str, abund_info: Optional[dict] = None) -> MyloContigID:
        abund_info = abund_info or {}
        header_string = header.split(" ")[0]
        match = cls.HEADER_PATTERN.match(header_string)

        if not match:
            raise ValueError(f"Unexpected contig header format for myloasm: '{header_string}'")

        groups = match.groupdict()
        return cls(
            _name=groups['name'],
            _abund_info=abund_info,
            length_val=int(groups['length']),
            is_circular=(groups['circular'] == 'yes'),
            mylo_depth=groups['depth'],
            is_duplicated=(groups['duplicated'] == 'yes')
        )

    @staticmethod
    def parse_name_from_header(header_string: str) -> str | None:
        """Parses just the contig name from a full header string."""
        match = MyloContigID.HEADER_PATTERN.match(header_string.split(" ")[0])
        if match:
            return match.group('name')
        return None

    @property
    def name(self) -> str: return self._name

    @property
    def length(self) -> int: return self.length_val


@dataclass
class GenericContigID(BaseContigID):
    """
    A minimal contig model where metadata is explicitly provided rather than inferred.
    """
    _provided_length: Optional[int] = None

    @classmethod
    def from_header(cls, header: str,
                    abund_info: Optional[dict] = None,
                    provided_length: Optional[int] = None) -> GenericContigID:
        abund_info = abund_info or {}
        raw = header.split(" ")[0]
        if raw.startswith(">"): raw = raw[1:]
        return cls(_name=raw, _abund_info=abund_info, _provided_length=provided_length)

    @staticmethod
    def parse_name_from_header(header_string: str) -> str | None:
        raw = header_string.split(" ")[0]
        if raw.startswith(">"): raw = raw[1:]
        return raw

    @property
    def name(self) -> str: return self._name

    @property
    def length(self) -> int | None: return self._provided_length


# Temporary aliases to maintain backward compatibility
ContigID = MyloContigID
_ParsedHeader = MyloParsedHeader

