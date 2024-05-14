import csv

def read_csv(csvfilename):
    """Reads a CSV file and returns the data as a tuple of tuples."""
    rows = ()
    with open(csvfilename) as csvfile:
        file_reader = csv.reader(csvfile)
        for row in file_reader:
            rows += (tuple(row), )
    return rows

def replicate(dna_strand):
    """Generates the complementary DNA strand by reversing and pairing bases."""
    dna_base_pairings = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }
    reversed_strand = dna_strand[::-1]
    ans = ""
    for base in reversed_strand:
        ans += dna_base_pairings[base]
    return ans

def find_transcription_region(dna_strand):
    """Locates and extracts the transcription region in a DNA strand using TATA and CGCG as markers."""
    tata_box = "TATA"
    stop_sequence = "CGCG"
    tata_index = dna_strand.find(tata_box)
    if tata_index == -1:
        return None
    start_index = tata_index + len(tata_box)
    stop_index = dna_strand.find(stop_sequence, start_index)
    if stop_index == -1:
        return None
    transcription_region = dna_strand[start_index:stop_index + len(stop_sequence)]
    return transcription_region

def transcribe(dna_strand):
    """Transcribes DNA to RNA within the transcription region."""
    dna_base_pairings = {
        "A": "U",
        "T": "A",
        "G": "C",
        "C": "G"
    }
    dna_strand = find_transcription_region(dna_strand)
    reversed_strand = dna_strand[::-1]
    ans = ""
    for base in reversed_strand:
        ans += dna_base_pairings[base]
    return ans

def reverse_transcribe(rna_strand):
    """Generates a complementary DNA strand from an RNA strand by reversing and pairing bases."""
    dna_base_pairings = {
        "A": "T",
        "U": "A",
        "G": "C",
        "C": "G"
    }
    reversed_strand = rna_strand[::-1]
    ans = ""
    for base in reversed_strand:
        ans += dna_base_pairings[base]
    return ans

def get_mapping(csvfilename):
    """Reads codon to amino acid mappings from a CSV file and returns a dictionary."""
    codon2aminomap = {}
    data = read_csv(csvfilename)
    for row in data:
        codon2aminomap[row[0]] = row[3]
    return codon2aminomap

def translate(rna_strand):
    """Translates an RNA strand to a protein sequence using codon mappings."""
    codon2amino = get_mapping("codon_mapping.csv")
    
    start_index = rna_strand.find('AUG')
    if start_index == -1:
        return None  

    protein = 'M'
    
    for i in range(start_index + 3, len(rna_strand), 3):
        codon = rna_strand[i:i+3]
        if len(codon) < 3:
            break  
        if codon in ['UAA', 'UAG', 'UGA']:
            protein += '_'
            return protein
        amino_acid = codon2amino.get(codon, '')
        protein += amino_acid

    return None
