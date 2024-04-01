
from flask import Flask, request, render_template

app = Flask(__name__)

def suffix_array(text):
    suffixes = [(text[i:], i) for i in range(len(text))]
    suffixes.sort(key=lambda x: x[0])
    return [suffix[1] for suffix in suffixes]

def bwt(text):
    return ''.join(text[i - 1] for i in suffix_array(text))

def create_fm_index(bwt):
    first_occurrence = {}
    count = {}
    for symbol in bwt:
        if symbol not in first_occurrence:
            first_occurrence[symbol] = 0
            count[symbol] = 0
    for i, symbol in enumerate(bwt):
        if i > 0:
            count[symbol] += 1
            if symbol not in count:
                count[symbol] = 0
    for symbol in sorted(count.keys()):
        if symbol not in first_occurrence:
            first_occurrence[symbol] = 0
    for symbol in sorted(count.keys()):
        if first_occurrence[symbol] == 0:
            first_occurrence[symbol] = count[symbol]
    return first_occurrence, count

def fm_substring_search(bwt, first_occurrence, count, pattern):
    occurrences = []
    top = 0
    bottom = len(bwt) - 1
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            if symbol in first_occurrence:
                top = first_occurrence[symbol] + count[symbol]
                if chr(ord(symbol) + 1) in first_occurrence:
                    bottom = first_occurrence[chr(ord(symbol) + 1)] - 1
            else:
                return occurrences
        else:
            for i in range(top, bottom + 1):
                occurrences.append(i)
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
    amino_acid_sequence = '-'.join(amino_acid_sequence)

    return protein_sequence, amino_acid_sequence



@app.route('/')
def index():
    return render_template('index.html')

@app.route('/suffix-array', methods=['POST'])
def suffix_array_search():
    text = request.form['text']
    suffixes = suffix_array(text)
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
    bwt_text = bwt(text)
    first_occurrence, count = create_fm_index(bwt_text)
    occurrences = fm_substring_search(bwt_text, first_occurrence, count, pattern)
    return render_template('index.html', occurrences=occurrences, text=text, pattern=pattern)


@app.route('/translate', methods=['POST'])
def translate():
    dna_sequence = request.form['dna_sequence']
    protein_sequence, amino_acid_sequence = translate_dna_to_protein(dna_sequence)
    return render_template('index.html', dna_sequence=dna_sequence, protein_sequence=protein_sequence, amino_acid_sequence=amino_acid_sequence)

if __name__ == '__main__':
    app.run(debug=True)

