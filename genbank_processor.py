#!/usr/bin/env python3
"""
GenBank Feature Table Processor - v4.0 (enhanced)
Converts GenBank files to NCBI-compliant feature table format (.tbl files)
with post-processing corrections for deprecated qualifiers and header formatting.

NEW in v4.0: Remove gb| prefix and | suffix from feature header
Enhanced: Robust sequence ID extraction, defensive cleaning, and edge-case handling
"""

import re
import os
import sys
from pathlib import Path


def extract_sequence_id(gb_content):
    """
    Extract a stable sequence ID for the feature table header.

    Priority:
      1. VERSION line (e.g., EU382073.1)
      2. ACCESSION line (first token if multiple)
      3. LOCUS line (second token)

    Returns a raw ID (without db prefix/suffix). Cleaning is applied separately.
    """
    version_re = re.compile(r'^VERSION\s+(\S+)', re.MULTILINE)
    accession_re = re.compile(r'^ACCESSION\s+(.+)$', re.MULTILINE)
    locus_re = re.compile(r'^LOCUS\s+(\S+)\s+(\S+)', re.MULTILINE)

    # Prefer VERSION
    m = version_re.search(gb_content)
    if m:
        version = m.group(1).strip()
        # Some records may have trailing tokens (rare). Keep only the first token.
        version = version.split()[0]
        return version

    # Fallback to ACCESSION
    m = accession_re.search(gb_content)
    if m:
        # ACCESSION can list multiple IDs; take the first token
        acc_line = m.group(1).strip()
        accession = acc_line.split()[0]
        return accession

    # Fallback to LOCUS name (second token is sequence length; we take first token as ID)
    m = locus_re.search(gb_content)
    if m:
        locus = m.group(1).strip()
        return locus

    # Last resort
    return "sequence"


def parse_features(gb_content):
    """
    Parse features from GenBank content.
    Extracts all feature types and their qualifiers.
    """
    features = []
    lines = gb_content.split('\n')
    in_features = False
    current_feature = None

    # Qualifiers to completely exclude from output
    excluded_qualifiers = {'codon_start', 'translation'}

    for i, line in enumerate(lines):
        if line.startswith('FEATURES'):
            in_features = True
            continue
        elif line.startswith('ORIGIN') or line.startswith('BASE COUNT'):
            in_features = False
            break

        if not in_features:
            continue

        # Skip the header line
        if 'Location/Qualifiers' in line:
            continue

        # Parse feature definition lines (5 spaces, then feature type)
        if line.startswith('     ') and not line.startswith('                   '):
            feature_match = re.match(r'     (\S+)\s+(.+)', line)
            if feature_match:
                feature_type = feature_match.group(1)
                location = feature_match.group(2).strip()

                if current_feature:
                    features.append(current_feature)

                current_feature = {
                    'type': feature_type,
                    'location': location,
                    'qualifiers': []
                }

        # Parse qualifier lines (19 spaces)
        elif line.startswith('                   ') and current_feature:
            # Accept word and dash in qualifier names (e.g., gene_synonym)
            qualifier_match = re.match(r'\s+/([\w\-]+)(?:=(.*))?', line)
            if qualifier_match:
                qual_name = qualifier_match.group(1)
                qual_value = qualifier_match.group(2)

                # Skip excluded qualifiers (codon_start, translation)
                if qual_name in excluded_qualifiers:
                    continue

                # Normalize and unwrap quoted values, support multi-line
                if qual_value is not None:
                    qual_value = qual_value.rstrip()

                    if qual_value.startswith('"'):
                        if qual_value.endswith('"') and len(qual_value) > 1:
                            qual_value = qual_value[1:-1]
                        else:
                            # Multi-line quoted value
                            qual_value = qual_value[1:]  # Remove opening quote
                            j = i + 1
                            while j < len(lines):
                                next_line = lines[j]
                                if next_line.startswith('                   ') and not next_line.strip().startswith('/'):
                                    continuation = next_line.strip()
                                    if continuation.endswith('"'):
                                        qual_value += " " + continuation[:-1]
                                        break
                                    else:
                                        qual_value += " " + continuation
                                    j += 1
                                else:
                                    break
                    else:
                        # Unquoted value: trim
                        qual_value = qual_value.strip()

                current_feature['qualifiers'].append({
                    'name': qual_name,
                    'value': qual_value if qual_value else None
                })

    # Add the last feature
    if current_feature:
        features.append(current_feature)

    return features


def _split_intervals(location):
    """
    Split a GenBank location into atomic intervals.
    Handles nested join/complement combinations.

    Returns list of interval strings like '123..456', '<1..>3311', or single '987'.
    """
    loc = location.strip()

    # Strip single-layer complement()
    is_complement = False
    if loc.startswith('complement(') and loc.endswith(')'):
        is_complement = True
        loc = loc[11:-1].strip()

    # If join(...) present, expand
    intervals = []
    if loc.startswith('join(') and loc.endswith(')'):
        inner = loc[5:-1]
        intervals = [p.strip() for p in inner.split(',')]
    else:
        intervals = [loc]

    return intervals, is_complement


def format_location_for_feature_table(location):
    """
    Convert GenBank location format to feature table format.
    Handles complement, join, and partial indicators (<, >).
    Also supports single-point intervals (e.g., '123').

    Returns list of 'start\\tend' strings, with partial markers preserved.
    """
    intervals, is_complement = _split_intervals(location)

    formatted_intervals = []
    for interval in intervals:
        interval = interval.strip()

        # Handle nested complement on each interval (e.g., join(complement(1..10),20..30))
        nested_complement = False
        if interval.startswith('complement(') and interval.endswith(')'):
            nested_complement = True
            interval = interval[11:-1].strip()

        # Single-point interval support
        if '..' in interval:
            start, end = interval.split('..', 1)
        else:
            start, end = interval, interval

        # Partial markers
        start_partial = '<' if start.startswith('<') else ''
        end_partial = '>' if end.startswith('>') else ''

        if start_partial:
            start = start[1:]
        if end_partial:
            end = end[1:]

        # Ensure numeric tokens remain strings; trust upstream parsing
        if (is_complement or nested_complement):
            # Reverse for complement intervals
            formatted_intervals.append(f"{start_partial}{end}\t{end_partial}{start}")
        else:
            formatted_intervals.append(f"{start_partial}{start}\t{end_partial}{end}")

    return formatted_intervals


def apply_qualifier_filters(features):
    """
    Apply post-processing filters to qualifiers:
    1. Remove 'gene' qualifier from mRNA features (already have gene feature)
    2. Remove 'label' qualifier from CDS features (deprecated qualifier)
    3. Change 'label' to 'note' in misc_feature (label is deprecated)
    """
    for feature in features:
        ftype = feature['type'].lower()

        # Rule 1
        if ftype == 'mrna':
            feature['qualifiers'] = [
                q for q in feature['qualifiers']
                if q['name'] != 'gene'
            ]

        # Rule 2
        if ftype == 'cds':
            feature['qualifiers'] = [
                q for q in feature['qualifiers']
                if q['name'] != 'label'
            ]

        # Rule 3
        elif ftype == 'misc_feature':
            for q in feature['qualifiers']:
                if q['name'] == 'label':
                    q['name'] = 'note'

    return features


def clean_sequence_id(sequence_id):
    """
    Remove database prefix and trailing pipe from sequence ID.

    Examples:
        gb|HMA4-1_Lan3.1_v2.1.0| -> HMA4-1_Lan3.1_v2.1.0
        gb|EU382073.1|           -> EU382073.1
        EU382073.1               -> EU382073.1
    """
    if sequence_id is None:
        return "sequence"

    seq = str(sequence_id).strip()

    # Remove a single leading "<db>|" prefix like gb|, ref|, emb|
    seq = re.sub(r'^[A-Za-z]+?\|', '', seq)

    # Remove trailing single pipe if present
    if seq.endswith('|'):
        seq = seq[:-1]

    # Final trim
    seq = seq.strip()

    # Guard against empty result
    if not seq:
        seq = "sequence"

    return seq


def write_feature_table(features, sequence_id, output_file):
    """
    Write features to feature table format (.tbl file).
    Output is NCBI-compliant with proper formatting.
    """
    clean_id = clean_sequence_id(sequence_id)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(f">Feature {clean_id}\n")

        for feature in features:
            # Skip source features
            if feature['type'].lower() == 'source':
                continue

            intervals = format_location_for_feature_table(feature['location'])

            # First interval with feature type
            f.write(f"{intervals[0]}\t{feature['type']}\n")

            # Additional intervals for multi-part features
            for interval in intervals[1:]:
                f.write(f"{interval}\n")

            # Qualifiers
            for qualifier in feature['qualifiers']:
                name = qualifier['name']
                value = qualifier['value']
                if value is not None and str(value).strip() != "":
                    f.write(f"\t\t\t{name}\t{value}\n")
                else:
                    f.write(f"\t\t\t{name}\n")


def process_genbank_file(input_file, output_dir=None):
    """
    Process a GenBank file and create feature table.

    Args:
        input_file: Path to input GenBank file
        output_dir: Directory to write output (default: same as input)

    Returns:
        Path to generated .tbl file
    """
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError:
        with open(input_file, 'r', encoding='latin-1') as f:
            content = f.read()

    # Extract sequence ID
    sequence_id = extract_sequence_id(content)

    # Parse features
    features = parse_features(content)

    # Apply post-processing filters
    features = apply_qualifier_filters(features)

    # Determine output file path
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    if output_dir:
        output_file = os.path.join(output_dir, f"{base_name}.tbl")
    else:
        output_file = f"{base_name}.tbl"

    # Write feature table
    write_feature_table(features, sequence_id, output_file)

    return output_file


def process_directory(input_dir, output_dir=None, pattern="*.gb"):
    """
    Process all GenBank files in a directory.

    Args:
        input_dir: Directory containing GenBank files
        output_dir: Directory for output files
        pattern: File pattern to match (default: *.gb)

    Returns:
        List of generated .tbl files
    """
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    input_path = Path(input_dir)
    output_files = []

    for gb_file in input_path.glob(pattern):
        if gb_file.is_file():
            try:
                output_file = process_genbank_file(str(gb_file), output_dir)
                output_files.append(output_file)
                print(f"✓ Processed: {gb_file.name} → {os.path.basename(output_file)}")
            except Exception as e:
                print(f"✗ Error processing {gb_file.name}: {e}")

    return output_files


def main():
    """Main entry point for command-line usage."""
    if len(sys.argv) < 2:
        print("GenBank Feature Table Processor v4.0 (enhanced)")
        print("\nUsage:")
        print("  python final_genbank_processor.py <input_file>")
        print("  python final_genbank_processor.py -i <input_dir> -o <output_dir>")
        print("  python final_genbank_processor.py -i <input_dir> -o <output_dir> --pattern '*.gbk'")
        print("\nExamples:")
        print("  # Process single file:")
        print("  python final_genbank_processor.py sequence.gb")
        print("\n  # Process directory with custom pattern:")
        print("  python final_genbank_processor.py -i ./genbank_files -o ./output")
        sys.exit(1)

    if sys.argv[1] in ['-i', '--input']:
        # Directory mode
        input_dir = sys.argv[2]
        output_dir = None
        pattern = "*.gb"

        if len(sys.argv) > 3 and sys.argv[3] in ['-o', '--output']:
            output_dir = sys.argv[4]

        if len(sys.argv) > 5 and sys.argv[5] == '--pattern':
            pattern = sys.argv[6]

        if not os.path.isdir(input_dir):
            print(f"Error: Input directory '{input_dir}' not found")
            sys.exit(1)

        output_files = process_directory(input_dir, output_dir, pattern)
        print(f"\nProcessing complete! Generated {len(output_files)} feature table files.")
    else:
        # Single file mode
        input_file = sys.argv[1]

        if not os.path.exists(input_file):
            print(f"Error: File '{input_file}' not found")
            sys.exit(1)

        try:
            output_file = process_genbank_file(input_file)
            print(f"✓ Successfully generated: {output_file}")
        except Exception as e:
            print(f"✗ Error processing file: {e}")
            sys.exit(1)


if __name__ == "__main__":
    main()
