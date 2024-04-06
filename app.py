from flask import Flask, request, render_template

app = Flask(__name__)

def build_suffix_array(text):
    """
    Build the suffix array for the given text.
    
    Args:
        text (str): The input text.
        
    Returns:
        list: The suffix array.
    """
    n = len(text)
    suffix_array = [i for i in range(n)]
    suffix_array.sort(key=lambda i: text[i:])
    return suffix_array

def bwt(text):
    """
    Compute the Burrows-Wheeler transform of the given text.
    
    Args:
        text (str): The input text.
        
    Returns:
        str: The Burrows-Wheeler transform of the text.
    """
    suffix_array = build_suffix_array(text)
    return ''.join(text[i - 1] for i in suffix_array)

def build_lcp_array(text, suffix_array):
    """
    Build the LCP (Longest Common Prefix) array for the given text and suffix array.
    
    Args:
        text (str): The input text.
        suffix_array (list): The suffix array.
        
    Returns:
        list: The LCP array.
    """
    n = len(text)
    lcp_array = [0] * n
    rank = [0] * n
    
    for i in range(n):
        rank[suffix_array[i]] = i
    
    l = 0
    for i in range(n):
        if rank[i] == n - 1:
            l = 0
            continue
        
        j = suffix_array[rank[i] + 1]
        while i + l < n and j + l < n and text[i + l] == text[j + l]:
            l += 1
        
        lcp_array[rank[i]] = l
        if l > 0:
            l -= 1
    
    return lcp_array

def find_first_occurrence(text, pattern, suffix_array, lcp_array):
    """
    Find the first occurrence of the pattern in the text using the FM-index.
    
    Args:
        text (str): The input text.
        pattern (str): The input pattern.
        suffix_array (list): The suffix array.
        lcp_array (list): The LCP array.
        
    Returns:
        int: The index of the first occurrence of the pattern in the text, or -1 if not found.
    """
    n = len(text)
    left, right = 0, n - 1
    
    while left <= right:
        mid = (left + right) // 2
        if lcp_array[mid] >= len(pattern):
            if text[suffix_array[mid]:suffix_array[mid] + len(pattern)] == pattern:
                return suffix_array[mid]
            elif text[suffix_array[mid]:suffix_array[mid] + len(pattern)] < pattern:
                left = mid + 1
            else:
                right = mid - 1
        elif text[suffix_array[mid]:suffix_array[mid] + lcp_array[mid]] < pattern[:lcp_array[mid]]:
            left = mid + 1
        else:
            right = mid - 1
    
    return -1

def count_occurrences(text, pattern, suffix_array, lcp_array):
    """
    Count the number of occurrences of the pattern in the text using the FM-index.
    
    Args:
        text (str): The input text.
        pattern (str): The input pattern.
        suffix_array (list): The suffix array.
        lcp_array (list): The LCP array.
        
    Returns:
        dict: A dictionary containing the first occurrence and count of each character in the pattern.
    """
    n = len(text)
    occurrences = {}
    
    for char in pattern:
        left, right = 0, n - 1
        first_occurrence = -1
        count = 0
        
        while left <= right:
            mid = (left + right) // 2
            if lcp_array[mid] >= len(char):
                if text[suffix_array[mid]:suffix_array[mid] + len(char)] == char:
              
                    left_idx = right_idx = mid
                    while left_idx > 0 and lcp_array[left_idx - 1] >= len(char):
                        left_idx -= 1
                    while right_idx < n - 1 and lcp_array[right_idx + 1] >= len(char):
                        right_idx += 1
                    first_occurrence = suffix_array[left_idx]
                    count = right_idx - left_idx + 1
                    break
                elif text[suffix_array[mid]:suffix_array[mid] + len(char)] < char:
                    left = mid + 1
                else:
                    right = mid - 1
            elif text[suffix_array[mid]:suffix_array[mid] + lcp_array[mid]] < char[:lcp_array[mid]]:
                left = mid + 1
            else:
                right = mid - 1
        
        occurrences[char] = (first_occurrence, count)
    
    return occurrences

def create_fm_index(text, pattern):
    """
    Create the FM-index for the given text and pattern.
    
    Args:
        text (str): The input text.
        pattern (str): The input pattern.
        
    Returns:
        dict: A dictionary containing the first occurrence and count of each character in the pattern.
    """
    suffix_array = build_suffix_array(text)
    lcp_array = build_lcp_array(text, suffix_array)

    occurrences = count_occurrences(text, pattern, suffix_array, lcp_array)
    
    return occurrences


def translate_dna_to_protein(dna_sequence):
    genetic_code = {
        'TTT': ('F', 'Phenylalanine'), 'TTC': ('F', 'Phenylalanine'),
        'TTA': ('L', 'Leucine'), 'TTG': ('L', 'Leucine'),
        'TCT': ('S', 'Serine'), 'TCC': ('S', 'Serine'), 'TCA': ('S', 'Serine'), 'TCG': ('S', 'Serine'),
        'TAT': ('Y', 'Tyrosine'), 'TAC': ('Y', 'Tyrosine'),
        'TGT': ('C', 'Cysteine'), 'TGC': ('C', 'Cysteine'), 'TGG': ('W', 'Tryptophan'),
        'CTT': ('L', 'Leucine'), 'CTC': ('L', 'Leucine'), 'CTA': ('L', 'Leucine'), 'CTG': ('L', 'Leucine'),
        'CCT': ('P', 'Proline'), 'CCC': ('P', 'Proline'), 'CCA': ('P', 'Proline'), 'CCG': ('P', 'Proline'),
        'CAT': ('H', 'Histidine'), 'CAC': ('H', 'Histidine'),
        'CAA': ('Q', 'Glutamine'), 'CAG': ('Q', 'Glutamine'),
        'CGT': ('R', 'Arginine'), 'CGC': ('R', 'Arginine'), 'CGA': ('R', 'Arginine'), 'CGG': ('R', 'Arginine'),
        'ATT': ('I', 'Isoleucine'), 'ATC': ('I', 'Isoleucine'), 'ATA': ('I', 'Isoleucine'),
        'ATG': ('M', 'Methionine'),
        'ACT': ('T', 'Threonine'), 'ACC': ('T', 'Threonine'), 'ACA': ('T', 'Threonine'), 'ACG': ('T', 'Threonine'),
        'AAT': ('N', 'Asparagine'), 'AAC': ('N', 'Asparagine'),
        'AAA': ('K', 'Lysine'), 'AAG': ('K', 'Lysine'),
        'AGT': ('S', 'Serine'), 'AGC': ('S', 'Serine'),
        'AGA': ('R', 'Arginine'), 'AGG': ('R', 'Arginine'),
        'GTT': ('V', 'Valine'), 'GTC': ('V', 'Valine'), 'GTA': ('V', 'Valine'), 'GTG': ('V', 'Valine'),
        'GCT': ('A', 'Alanine'), 'GCC': ('A', 'Alanine'), 'GCA': ('A', 'Alanine'), 'GCG': ('A', 'Alanine'),
        'GAT': ('D', 'Aspartic Acid'), 'GAC': ('D', 'Aspartic Acid'),
        'GAA': ('E', 'Glutamic Acid'), 'GAG': ('E', 'Glutamic Acid'),
        'GGT': ('G', 'Glycine'), 'GGC': ('G', 'Glycine'), 'GGA': ('G', 'Glycine'), 'GGG': ('G', 'Glycine'),
        'TAA': ('*', 'Stop'), 'TAG': ('*', 'Stop'), 'TGA': ('*', 'Stop')
    }

    if dna_sequence is not None:
        dna_sequence = dna_sequence.upper().replace(' ', '')
        protein_sequence = []
        amino_acid_sequence = []

        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i + 3]
            if len(codon) < 3:
                break  # Skip incomplete codons at the end of the sequence

            if codon in genetic_code:
                amino_acid, aa_name = genetic_code[codon]
                protein_sequence.append(amino_acid)
                amino_acid_sequence.append(aa_name)
            else:
                amino_acid_sequence.append('Unknown')

        protein_sequence = ''.join(protein_sequence)
        amino_acid_sequence = ' - '.join(amino_acid_sequence)

        return protein_sequence, amino_acid_sequence



    if dna_sequence is not None:
        dna_sequence = dna_sequence.upper().replace(' ', '')
        protein_sequence = []
        amino_acid_sequence = []

        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i + 3]
            if len(codon) < 3:
                break  # Skip incomplete codons at the end of the sequence

            if codon in genetic_code:
                amino_acid, aa_name = genetic_code[codon]
                protein_sequence.append(amino_acid)
                amino_acid_sequence.append(aa_name)
            else:
                amino_acid_sequence.append('Unknown')

        protein_sequence = ''.join(protein_sequence)
        amino_acid_sequence = ' - '.join(amino_acid_sequence)

        return protein_sequence, amino_acid_sequence

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/suffix-array', methods=['POST'])
def suffix_array_search():
    text = request.form['text']
    suffixes = build_suffix_array(text)
    return render_template('index.html', suffixes=suffixes, text=text)

@app.route('/bwt', methods=['POST'])
def bwt_search():
    text = request.form['text']
    bwt_text = bwt(text)
    return render_template('index.html', bwt_text=bwt_text, text=text)

@app.route('/fm-index', methods=['POST'])
def fm_index():
    text = request.form['text']
    pattern = request.form['pattern']
    occurrences = create_fm_index(text, pattern)


    result = "FM Index:\n"
    result += "{\n"
    for char, (first_occurrence, count) in occurrences.items():
        result += f"    '{char}': {{'First Occurrence': {first_occurrence}, 'Count': {count}}},\n"
    result += "}\n"

    return render_template('index.html', result=result)



@app.route('/dna-to-protein', methods=['POST'])
def dna_to_protein():
    dna_sequence = request.form['dna_sequence']
    protein_sequence, amino_acid_sequence = translate_dna_to_protein(dna_sequence)
    return render_template('index.html', protein_sequence=protein_sequence, amino_acid_sequence=amino_acid_sequence, dna_sequence=dna_sequence)

if __name__ == '__main__':
    app.run(debug=True)
