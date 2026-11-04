"""
Problem: Number System Conversion

Implement conversions between binary, decimal and hexadecimal,
as well as a general function for converting between arbitrary
number systems (bases 2 to 36).

Functions:
- binary_to_decimal()
- binary_to_hexadecimal()
- decimal_to_hexadecimal()
- hexadecimal_to_decimal()
- hexadecimal_to_binary()
- decimal_to_binary()
- convert_base(): General conversion between arbitrary bases.

Approach:
Two helper functions are used:

1. _to_decimal():
   Converts a number in an arbitrary base into its decimal
   representation using positional notation.

2. _from_decimal():
   Converts a decimal integer into an arbitrary base using
   repeated division and remainders.

The general convert_base() function combines these two
operations.

Digits 0-9 and A-Z are supported. Output uses uppercase
letters and removes unnecessary leading zeros.

Complexity Analysis:
Let d denote the number of digits in the input or output.

- Time Complexity: O(d^2) bit operations.
  The conversion uses repeated arithmetic operations on
  arbitrary-precision integers. Each operation can take
  O(d) time, and there are O(d) iterations.

  Under a unit-cost arithmetic model, the complexity is O(d).

- Space Complexity: O(d).
  The intermediate integer and output representation require
  space proportional to the number of digits.

Python Alternatives:
- int(binary, 2): Converts binary to decimal.
- bin(decimal): Converts decimal to binary.
- hex(decimal): Converts decimal to hexadecimal.
- int(hexadecimal, 16): Converts hexadecimal to decimal.
- format(decimal, 'b'): Converts decimal to binary.
- format(decimal, 'X'): Converts decimal to uppercase hex.
"""


_DIGITS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _validate_base(base: int) -> None:
    if type(base) is not int:
        raise TypeError("Base must be an integer.")

    if not 2 <= base <= 36:
        raise ValueError("Base must be between 2 and 36.")


def _to_decimal(number: str, base: int) -> int:
    """Convert a number from an arbitrary base to decimal."""

    if not isinstance(number, str):
        raise TypeError("Number must be a string.")

    number = number.strip().upper()

    if not number:
        raise ValueError("Number cannot be empty.")

    sign = 1

    if number[0] == "-":
        sign = -1
        number = number[1:]

    elif number[0] == "+":
        number = number[1:]

    if not number:
        raise ValueError("Invalid number.")

    result = 0

    for character in number:
        digit = _DIGITS.find(character)

        if digit == -1 or digit >= base:
            raise ValueError(f"Invalid digit '{character}' for base {base}.")

        result = result * base + digit

    return sign * result


def _from_decimal(number: int, base: int) -> str:
    """Convert a decimal integer to an arbitrary base."""

    if type(number) is not int:
        raise TypeError("Number must be an integer.")

    if number == 0:
        return "0"

    sign = ""

    if number < 0:
        sign = "-"
        number = -number

    result = []

    while number > 0:
        number, remainder = divmod(number, base)
        result.append(_DIGITS[remainder])

    result.reverse()

    return sign + "".join(result)


# --------------------------------------------------
# General conversion between arbitrary number systems
# --------------------------------------------------

def convert_base(number: str, from_base: int, to_base: int) -> str:
    """
    Convert a number from one arbitrary base to another.

    Example:
    convert_base("1F", 16, 2) returns "11111".
    """

    _validate_base(from_base)
    _validate_base(to_base)

    decimal = _to_decimal(number, from_base)

    return _from_decimal(decimal, to_base)


# --------------------------------------------------
# Binary conversions
# --------------------------------------------------

def binary_to_decimal(binary: str) -> int:
    return _to_decimal(binary, 2)


def binary_to_hexadecimal(binary: str) -> str:
    return convert_base(binary, 2, 16)


# --------------------------------------------------
# Decimal conversions
# --------------------------------------------------

def decimal_to_hexadecimal(decimal: int) -> str:
    return _from_decimal(decimal, 16)


def decimal_to_binary(decimal: int) -> str:
    return _from_decimal(decimal, 2)


# --------------------------------------------------
# Hexadecimal conversions
# --------------------------------------------------

def hexadecimal_to_decimal(hexadecimal: str) -> int:
    return _to_decimal(hexadecimal, 16)


def hexadecimal_to_binary(hexadecimal: str) -> str:
    return convert_base(hexadecimal, 16, 2)
