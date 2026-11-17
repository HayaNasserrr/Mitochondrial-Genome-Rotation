"""Rotate an annotated circular mitochondrial genome to the start of tRNA-Pro.

Requires a FASTA file and a matching annotated GenBank file.
"""
from copy import deepcopy
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SimpleLocation, CompoundLocation


def rotate_location(location, cut, genome_length):
    """Shift a Biopython location, retaining biological part order and strand.

    Split each interval at the new origin, then merge adjacent pieces. Works
    with simple locations, compound locations, and either strand orientation.
    """
    if location is None:
        raise ValueError("An annotation has no location.")
    if isinstance(location, CompoundLocation) and location.operator != "join":
        raise ValueError(f"Unsupported compound location operator: {location.operator}")

    shifted_parts = []
    for part in location.parts:
        start, end = int(part.start), int(part.end)
        strand = part.strand
        if not 0 <= start < end <= genome_length:
            raise ValueError(f"Invalid feature interval: {part}")

        # Two possible segments: before the new origin and after it.
        # Each segment is shifted independently, without extending across
        # the entire circular genome.
        segments = []
        if start < cut:
            left_end = min(end, cut)
            if start < left_end:
                segments.append(SimpleLocation(start + genome_length - cut,
                                               left_end + genome_length - cut,
                                               strand=strand))
        if end > cut:
            right_start = max(start, cut)
            if right_start < end:
                segments.append(SimpleLocation(right_start - cut,
                                               end - cut, strand=strand))

        # A negative-strand feature is extracted in reverse part order.
        # Preserve the nucleotide sequence, not merely the coordinate set.
        if strand == -1:
            segments.reverse()  # reverse biological order for - strand
        shifted_parts.extend(segments)

    # Merge neighboring intervals in biological extraction order.
    merged = []
    for part in shifted_parts:
        if merged and merged[-1].strand == part.strand:
            previous = merged[-1]
            if part.strand == -1 and int(part.end) == int(previous.start):
                merged[-1] = SimpleLocation(part.start, previous.end,
                                            strand=part.strand)
                continue
            if part.strand != -1 and int(previous.end) == int(part.start):
                merged[-1] = SimpleLocation(previous.start, part.end,
                                            strand=part.strand)
                continue
        merged.append(part)

    return merged[0] if len(merged) == 1 else CompoundLocation(merged, operator="join")


def run():
    print("Mitochondrial Genome Rotation Tool")
    fasta_path = input("Enter FASTA file name: ").strip().strip('"')
    genbank_path = input("Enter GenBank file name: ").strip().strip('"')
    fasta_record = SeqIO.read(fasta_path, "fasta")
    original_record = SeqIO.read(genbank_path, "genbank")
    original_seq = fasta_record.seq
    genome_length = len(original_seq)

    print("\nInput validation:")
    if genome_length == 0 or original_seq.upper() != original_record.seq.upper():
        raise ValueError("FASTA and GenBank sequences do not match, or are empty.")
    print("Sequence validation: PASSED")
    print("Genome length:", genome_length)

    trnp_features = [f for f in original_record.features
                     if f.type == "gene" and
                     f.qualifiers.get("gene", [""])[0].upper() in {"TRNP", "TRNP(TGG)", "TRNP(UGG)"}]
    if len(trnp_features) != 1:
        raise ValueError("Exactly one annotated TRNP gene is required.")
    trnp = trnp_features[0]
    if trnp.location is None:
        raise ValueError("TRNP has no location.")

    # On a circular molecule, the 5'-most coordinate of the annotation is
    # represented by the first interval in biological part order. For the
    # project specification, rotation starts at its lowest coordinate when
    # the feature is contiguous. Wrapped features require an explicit origin.
    if isinstance(trnp.location, CompoundLocation):
        raise ValueError("TRNP spans the current origin; its start is ambiguous."
                         " Please supply a GenBank annotation with an unwrapped TRNP gene.")
    cut = int(trnp.location.start)
    print("\nTRNP found: PASSED")
    print("TRNP start:", cut + 1)
    print("TRNP end:", int(trnp.location.end))
    print("TRNP strand:", trnp.location.strand)

    if cut == 0:
        print("\nThe mitochondrial genome already starts at tRNAPro.")
        print("No rotation is needed.")
        return
    print("\nThe mitochondrial genome does not start at tRNAPro.")
    print("Rotation is needed.")
    while True:
        answer = input("Would you like to rotate the sequence to start at tRNAPro? (Y/N): ").strip().upper()
        if answer in {"Y", "N"}:
            break
        print("Invalid choice. Please enter Y or N.")
    if answer == "N":
        print("Rotation cancelled. Original sequence unchanged.")
        return

    rotated_seq = original_seq[cut:] + original_seq[:cut]
    print("\nRotation completed.")
    print("Original length:", genome_length)
    print("Rotated length:", len(rotated_seq))
    print("Rotation index:", cut)
    if len(rotated_seq) != genome_length or rotated_seq[genome_length-cut:] + rotated_seq[:genome_length-cut] != original_seq:
        raise ValueError("Circular sequence validation failed.")
    print("Circular sequence validation: PASSED")

    rotated_record = deepcopy(original_record)
    rotated_record.seq = rotated_seq
    print("\nUpdating GenBank annotations:")
    for original_feature, updated_feature in zip(original_record.features, rotated_record.features):
        if original_feature.type == "source":
            # Source spans the entire circular molecule, regardless of origin.
            updated_feature.location = SimpleLocation(0, genome_length,
                                                       strand=original_feature.location.strand)
        else:
            updated_feature.location = rotate_location(original_feature.location, cut, genome_length)

    print("\nValidating all GenBank features:")
    validated = 0
    for original_feature, updated_feature in zip(original_record.features, rotated_record.features):
        if original_feature.type == "source":
            continue
        label = original_feature.qualifiers.get("gene", ["Unknown"])[0]
        context = f"{original_feature.type} ({label})"
        if len(original_feature) != len(updated_feature):
            raise ValueError(f"Feature length changed: {context}; old={original_feature.location}, new={updated_feature.location}")
        if original_feature.location.strand != updated_feature.location.strand:
            raise ValueError(f"Feature strand changed: {context}")
        if original_feature.extract(original_seq).upper() != updated_feature.extract(rotated_seq).upper():
            raise ValueError(f"Feature DNA sequence changed: {context}")
        validated += 1
    print("Validated features:", validated)
    print("Feature validation: PASSED")

    trnp_after = [f for f in rotated_record.features if f.type == "gene"
                  and f.qualifiers.get("gene", [""])[0].upper() in {"TRNP", "TRNP(TGG)", "TRNP(UGG)"}]
    if len(trnp_after) != 1 or int(trnp_after[0].location.start) != 0:
        raise ValueError("TRNP did not move to position 1.")

    rotated_fasta = deepcopy(fasta_record)
    rotated_fasta.seq = rotated_seq
    rotated_fasta.id = fasta_record.id + "_rotated"
    rotated_fasta.description = rotated_fasta.id
    rotated_record.id = original_record.id + "_rotated"
    rotated_record.name = (original_record.name or "mtDNA")[:12] + "_rot"
    rotated_record.description = "Mitochondrial genome rotated to TRNP start"
    rotated_record.annotations["topology"] = "circular"
    rotated_record.annotations.pop("accessions", None)
    rotated_record.annotations.pop("sequence_version", None)

    output_fasta = Path("rotated_genome.fasta")
    output_gb = Path("rotated_genome.gb")
    SeqIO.write(rotated_fasta, output_fasta, "fasta")
    SeqIO.write(rotated_record, output_gb, "genbank")

    saved_fasta = SeqIO.read(output_fasta, "fasta")
    saved_gb = SeqIO.read(output_gb, "genbank")
    if saved_fasta.seq != rotated_seq or saved_gb.seq != rotated_seq:
        raise ValueError("Saved sequence validation failed.")
    if len(saved_gb.features) != len(original_record.features):
        raise ValueError("Saved feature count changed.")
    for before, after in zip(rotated_record.features, saved_gb.features):
        if before.type != after.type or len(before) != len(after) or before.extract(rotated_seq).upper() != after.extract(rotated_seq).upper():
            raise ValueError(f"GenBank round-trip validation failed: {before.type}")
    print("\nSaved files validation: PASSED")
    print("Saved TRNP start: 1")
    print("Saved feature count:", len(saved_gb.features))
    print("FASTA saved:", output_fasta)
    print("GenBank saved:", output_gb)
    print("Files saved in:", Path.cwd())


if __name__ == "__main__":
    run()
