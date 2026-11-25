"""Maestro CT reader adapted from feynman-ai's tool_server/app/tools/mae.py."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
import re


@dataclass
class MaeTable:
    columns: list[str]
    rows: dict[int, dict[str, str | None]]


@dataclass
class MaeBlock:
    name: str
    properties: dict[str, str | None]
    tables: dict[str, MaeTable]


class MAEParser:
    _separator = ":::"
    _empty_value = "<>"
    _component = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\[(\d+)\]")
    _property = re.compile(r"[birs]_[A-Za-z0-9_]+")

    class Tokenizer:
        _pattern = re.compile(
            r'"(?:[^"\\]|\\.)*"|#[^\n]*|:::|[{}]|[^\s"{}#]+|\S', re.DOTALL
        )

        def __init__(self, text):
            self._matches = (
                match.group() for match in self._pattern.finditer(text)
                if not match.group().startswith("#")
            )
            self.current = next(self._matches, None)

        def read_token(self):
            current = self.current
            self.current = next(self._matches, None)
            return current

    def parse(self, path):
        path = Path(path)
        opener = gzip.open if path.suffix.lower() in {".maegz", ".gz"} else open
        with opener(path, "rt", encoding="utf-8-sig") as source:
            self.tokenizer = self.Tokenizer(source.read())

        entries = []
        while self.tokenizer.current is not None:
            entries.append(self.parse_entry(self.accept_block()))
        return entries

    def accept_block(self):
        if self.tokenizer.current == "{":
            self.expect("{")
            return ""
        name = self.tokenizer.read_token()
        if name not in {"f_m_ct", "p_m_ct"}:
            raise ValueError(f"Unexpected MAE block {name!r}")
        self.expect("{")
        return name

    def parse_entry(self, name):
        properties = self.parse_properties()
        components = {}
        while self.tokenizer.current != "}":
            component = self.accept_component()
            if component is None:
                continue
            kind, count = component
            if kind in components:
                raise ValueError(f"Duplicate MAE table {kind}")
            components[kind] = self.parse_component(kind, count)
        self.expect("}")
        return MaeBlock(name, properties, components)

    def accept_component(self):
        token = self.tokenizer.read_token()
        if token is None:
            raise ValueError("Unexpected end of MAE block")
        match = self._component.fullmatch(token)
        self.expect("{")
        if match is None:
            self.ignore_block()
            return None
        return match[1], int(match[2])

    def ignore_block(self):
        depth = 1
        while depth:
            token = self.tokenizer.read_token()
            if token is None:
                raise ValueError("Unterminated MAE subblock")
            if token == "{":
                depth += 1
            elif token == "}":
                depth -= 1

    def parse_properties(self):
        columns = self.parse_component_columns()
        return self.read_values(columns)

    def parse_component_columns(self):
        columns = []
        while self.tokenizer.current != self._separator:
            name = self.tokenizer.read_token()
            if name is None or not self._property.fullmatch(name):
                raise ValueError(f"Invalid MAE property name {name!r}")
            columns.append(name)
        self.expect(self._separator)
        if len(columns) != len(set(columns)):
            raise ValueError("Duplicate MAE property name")
        return columns

    def parse_component(self, name, count):
        columns = self.parse_component_columns()
        rows = {}
        for _ in range(count):
            token = self.tokenizer.read_token()
            if token is None or not token.isdecimal():
                raise ValueError(f"Invalid {name} row index {token!r}")
            index = int(token)
            if index < 1 or index in rows:
                raise ValueError(f"Invalid or duplicate {name} row index {index}")
            rows[index] = self.read_values(columns)
        self.expect(self._separator)
        self.expect("}")
        return MaeTable(columns, rows)

    def read_values(self, columns):
        row = {}
        for name in columns:
            token = self.tokenizer.read_token()
            if token in {None, "{", "}", self._separator, '"'}:
                raise ValueError(f"Missing or malformed MAE value for {name}")
            if token.startswith('"'):
                row[name] = re.sub(r'\\(["\\])', r'\1', token[1:-1])
            else:
                row[name] = None if token == self._empty_value else token
        return row

    def expect(self, expected):
        token = self.tokenizer.read_token()
        if token != expected:
            raise ValueError(f"Expected {expected!r} in MAE, got {token!r}")
