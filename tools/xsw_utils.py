"""Helpers for reading and updating parameters in TransportCR XSW files."""

from __future__ import print_function

import os
import re


class XswParameterError(ValueError):
    """Raised when an XSW parameter cannot be read unambiguously."""


def _parameter_tag_pattern(parameter_name):
    escaped_name = re.escape(parameter_name)
    return re.compile(
        r'<[^!?/][^>]*\btitle\s*=\s*["\']'
        + escaped_name
        + r'["\'][^>]*>',
        re.DOTALL,
    )


_VALUE_PATTERN = re.compile(r'\bvalue\s*=\s*(["\'])(.*?)\1', re.DOTALL)
_INVALID_PARAMETER_VALUE_CHARACTERS = frozenset('&<"\'')


def validate_parameter_value(value):
    """Reject characters that cannot be inserted as a plain XML attribute."""
    value = str(value)
    invalid = sorted(set(value) & _INVALID_PARAMETER_VALUE_CHARACTERS)
    if invalid:
        rendered = ", ".join(repr(character) for character in invalid)
        raise XswParameterError(
            "parameter value contains unsupported XML character(s): {0}".format(
                rendered
            )
        )
    return value


def parameter_values_from_text(text, parameter_name):
    """Return every value attribute for tags with the requested title."""
    values = []
    for tag_match in _parameter_tag_pattern(parameter_name).finditer(text):
        value_match = _VALUE_PATTERN.search(tag_match.group(0))
        if value_match:
            values.append(value_match.group(2))
    return values


def get_parameter_values(file_name, parameter_name):
    with open(file_name, "r") as xsw_file:
        return parameter_values_from_text(xsw_file.read(), parameter_name)


def get_parameter(file_name, parameter_name):
    """Read exactly one parameter value or raise XswParameterError."""
    values = get_parameter_values(file_name, parameter_name)
    if not values:
        raise XswParameterError(
            "parameter {0} not found in {1}".format(parameter_name, file_name)
        )
    if len(values) != 1:
        raise XswParameterError(
            "parameter {0} occurs {1} times in {2}".format(
                parameter_name, len(values), file_name
            )
        )
    return values[0]


def replace_parameter_in_text(text, parameter_name, new_value):
    """Replace value attributes and return (updated_text, replacement_count)."""
    new_value = validate_parameter_value(new_value)
    replacement_count = 0

    def replace_tag(tag_match):
        nonlocal replacement_count
        tag = tag_match.group(0)

        def replace_value(value_match):
            nonlocal replacement_count
            replacement_count += 1
            quote = value_match.group(1)
            return "value={0}{1}{0}".format(quote, new_value)

        return _VALUE_PATTERN.sub(replace_value, tag, count=1)

    updated = _parameter_tag_pattern(parameter_name).sub(replace_tag, text)
    return updated, replacement_count


def write_text_atomic(file_name, text):
    """Write text through the historical .tmp sibling and replace atomically."""
    tmp_file_name = str(file_name) + ".tmp"
    with open(tmp_file_name, "w") as output_file:
        output_file.write(text)
    os.replace(tmp_file_name, file_name)


def replace_parameter(file_name, parameter_name, new_value):
    """Replace exactly one parameter in a file and return the match count."""
    with open(file_name, "r") as xsw_file:
        text = xsw_file.read()
    updated, replacement_count = replace_parameter_in_text(
        text, parameter_name, str(new_value)
    )
    if replacement_count == 1:
        write_text_atomic(file_name, updated)
    return replacement_count
