def is_in_gc_bounds(seq: str, gc_bounds: int | float | tuple = (0, 100)) -> bool:
    seq = seq.lower()
    gc_count = 100 * (seq.count('g') + seq.count('c')) / len(seq)

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    return gc_bounds[0] <= gc_count <= gc_bounds[1]


def is_in_length_bounds(seq: str, length_bounds: int | float | tuple = (0, 2**32)) -> bool:
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def is_above_quality_threshold(quality: str, quality_threshold: float = 0) -> bool:
    quality_mean = sum([ord(letter) - 33 for letter in quality]) / len(quality)
    return quality_mean >= quality_threshold


def is_read_good(seq: str, quality: str,
                 gc_bounds: int | float | tuple = (0, 100),
                 length_bounds: int | float | tuple = (0, 2 ** 32),
                 quality_threshold: int | float = 0) -> bool:
    if (
            not is_in_gc_bounds(seq, gc_bounds) or
            not is_in_length_bounds(seq, length_bounds) or
            not is_above_quality_threshold(quality, quality_threshold)
    ):
        return False

    return True