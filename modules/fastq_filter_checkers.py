def is_in_gc_bounds(seq: str, gc_bounds: int | float | tuple = (0, 100)) -> bool:
    seq = seq.lower()
    gc_count = 100 * (seq.count('g') + seq.count('c')) / len(seq)

    if type(gc_bounds) is not tuple:
        gc_bounds = (0, gc_bounds)

    if gc_bounds[0] <= gc_count <= gc_bounds[1]:
        return True

    return False


def is_in_length_bounds(seq: str, length_bounds: int | float | tuple = (0, 2**32)) -> bool:
    if type(length_bounds) is not tuple:
        length_bounds = (0, length_bounds)

    if length_bounds[0] <= len(seq) <= length_bounds[1]:
        return True

    return False


def is_above_quality_threshold(quality: str, quality_threshold: float = 0) -> bool:
    quality_sum=0

    for l in quality:
        quality_sum += ord(l)-33
    quality_mean = quality_sum / len(quality)

    if quality_mean >= quality_threshold:
        return True

    return False


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