class BaseMaskConfigException(Exception):
    """Raise if base mask for demux_reads is malformed"""

    pass


def split_bases_mask(bases_mask):
    """Parse a picard style bases mask and return a list of tuples
    >>> split_bases_mask("100T8B8B100T")
    [('T', 100), ('B', 8), ('B', 8), ('T', 100)]
    """
    splitat = []
    for i, c in enumerate(bases_mask):
        if c.isalpha():
            splitat.append(i)

    # Check that mask is well-behaved
    if splitat[0] == 0:
        raise BaseMaskConfigException("Mask must start with number of cycles, not type")
    # Check that no letters appear next to each other
    diffs = []
    for i in range(len(splitat) - 1):
        diffs.append(splitat[i + 1] - splitat[i])
    if (0 in diffs) or (1 in diffs):
        raise BaseMaskConfigException("Type characters must be separated by a number (of cycles)")

    result = []
    num = ""
    for i, _character in enumerate(bases_mask):
        if i not in splitat:
            num += bases_mask[i]
        elif int(num) == 0:
            pass
        else:
            result.append((bases_mask[i].upper(), int(num)))
            num = ""

    return result


def compare_bases_mask(planned_reads, bases_mask, demux_tool="bcl2fastq"):
    """Match user input bases mask to planned_reads from flowcell and decide if compatible.
    Return list of lists of tuples with type and number of cycles"""

    planned = split_bases_mask(planned_reads)
    mask = split_bases_mask(bases_mask)

    lengths1 = [count for _type, count in planned]
    lengths2 = [count for _type, count in mask]
    if not sum(lengths1) == sum(lengths2):
        raise BaseMaskConfigException("Your base mask has more or fewer cycles than planned")

    if demux_tool != "bcl2fastq":
        return [[x] for x in mask]

    exp_mask = []
    for op, cycles in mask:
        for i in range(cycles):
            exp_mask.append((op, 1))

    result = []
    offset = 0
    for p_type, p_cycles in planned:
        curr = []
        for m_type, _ in exp_mask[offset : offset + p_cycles]:
            if curr and curr[-1][0] == m_type:
                curr[-1][1] += 1
            else:
                curr.append([m_type, 1])
        result.append(list(map(tuple, curr)))
        offset += p_cycles

    return result


def translate_tuple_to_basemask(tup, demux_tool):
    """Return illumina or picard-style base mask string"""

    picard_to_illumina = {"T": "y", "S": "n", "B": "I"}

    if demux_tool == "bcl2fastq":
        return picard_to_illumina[tup[0]] + str(tup[1])
    else:
        return str(tup[1]) + tup[0]


def return_bases_mask(planned_reads, demux_reads, demux_tool="bcl2fastq"):
    """Parse planned_reads and demux_reads (user-configured base mask), compare for compatiblity
    and return a string to either give to bcl2fastq or picard"""

    if "M" in demux_reads and demux_tool == "bcl2fastq":
        raise BaseMaskConfigException("You cannot assign UMIs ('M') if using bcl2fastq")

    mask_list = compare_bases_mask(planned_reads, demux_reads, demux_tool)

    new_mask = []
    for lst in mask_list:
        substr = [translate_tuple_to_basemask(t, demux_tool) for t in lst]
        new_mask.append("".join(substr))
    return ",".join(new_mask) if demux_tool == "bcl2fastq" else "".join(new_mask)
