"""The centesimal division of Bouvard's 1821 tables, converted to the sexagesimal one.

Bouvard (p. xviii): the circle has 400 degrees, the degree 100 minutes, the minute 100
seconds; the day has 10 hours, the hour 100 minutes, the minute 100 seconds. A centesimal
degree is 0.9 of a sexagesimal degree, and a centesimal second is 0.324 of an arcsecond:
0.9 * 3600 / 10000.
"""

from __future__ import annotations

GRADE_IN_DEG = 0.9  # 360 / 400
CENTESIMAL_SEC_IN_ARCSEC = 0.324  # 360 * 3600 / (400 * 100 * 100)
DECIMAL_HOUR_IN_DAY = 0.1


def grades_to_deg(grades: float) -> float:
    """Centesimal degrees (400 to the circle) to sexagesimal degrees."""
    return grades * GRADE_IN_DEG


def centesimal_dms_to_deg(degrees: int, minutes: int, seconds: float) -> float:
    """A centesimal angle printed as degrees, minutes (0 to 99) and seconds (below 100), to
    sexagesimal degrees. The sign of ``degrees`` is not applied to the other parts: the
    tables print only positive angles."""
    if degrees < 0 or not 0 <= minutes < 100 or not 0 <= seconds < 100:
        raise ValueError(f"not a centesimal angle: {degrees} {minutes} {seconds}")
    return grades_to_deg(degrees + minutes / 100 + seconds / 10000)


def centesimal_sec_to_arcsec(seconds: float) -> float:
    """Centesimal seconds (1,000,000 to the right angle) to arcseconds."""
    return seconds * CENTESIMAL_SEC_IN_ARCSEC


def decimal_time_to_day(hours: int, minutes: int, seconds: float) -> float:
    """Decimal time of day (10 hours of 100 minutes of 100 seconds) to a fraction of a day."""
    if not 0 <= hours < 10 or not 0 <= minutes < 100 or not 0 <= seconds < 100:
        raise ValueError(f"not a decimal time: {hours} {minutes} {seconds}")
    return (hours + minutes / 100 + seconds / 10000) * DECIMAL_HOUR_IN_DAY
