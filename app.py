from flask import Flask, render_template, request, jsonify
import io
import os
from io import StringIO
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
import pandas as pd
import numpy as np


app = Flask(__name__)


def perform_global_alignment(seq1, seq2, match_score, mismatch_score, open_gap_score, extend_gap_score):
    if not seq1 or not seq2:
        return "Error: Empty input sequences."

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = match_score      
    aligner.mismatch_score = mismatch_score 
    aligner.open_gap_score = open_gap_score  
    aligner.extend_gap_score = extend_gap_score 

    try:
        alignments = aligner.align(seq1, seq2)
        if not alignments:
            return "No alignment possible."
        best_alignment = alignments[0]
        aligned_seq1, aligned_seq2 = best_alignment[0], best_alignment[1]
        score = best_alignment.score
        
        symbols = []
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == b:
                symbols.append("|")
            else:
                symbols.append(".")
        
        
        formatted_output = (
            f"Aligned Sequence 1: {aligned_seq1}\n"
            f"                   {''.join(symbols)}\n"
            f"Aligned Sequence 2: {aligned_seq2}\n"
            f"Alignment Score: {score}"
        )
        return formatted_output

    except Exception as e:
        return f"Alignment error: {str(e)}"

def build_suffix_array(text):
    n = len(text)
    suffix_array = [i for i in range(n)]
    suffix_array.sort(key=lambda i: text[i:])
    return suffix_array

def search_with_suffix_array(text, pattern):
    suffix_array = build_suffix_array(text)
    n = len(text)
    m = len(pattern)
    occurrences = []

    left, right = 0, n - 1
    while left <= right:
        mid = (left + right) // 2
        suffix = text[suffix_array[mid]:suffix_array[mid] + m]
        if suffix == pattern:
            start = mid
            while start >= 0 and text[suffix_array[start]:suffix_array[start] + m] == pattern:
                occurrences.append(suffix_array[start])
                start -= 1
            end = mid + 1
            while end < n and text[suffix_array[end]:suffix_array[end] + m] == pattern:
                occurrences.append(suffix_array[end])
                end += 1
            break
        elif suffix < pattern:
            left = mid + 1
        else:
            right = mid - 1

    return sorted(occurrences)

def bwt(text):
    suffix_array = build_suffix_array(text)
    return ''.join(text[i - 1] for i in suffix_array)

def create_fm_index(text, pattern):
    bwt_text = bwt(text)
    sorted_bwt = sorted(bwt_text)
    occurrences = {}
    first_occurrence = {}
    current_char = None
    for idx, char in enumerate(sorted_bwt):
        if char != current_char:
            first_occurrence[char] = idx
            current_char = char
    unique_chars_in_pattern = set(pattern)
    for char in unique_chars_in_pattern:
        fo = first_occurrence.get(char, -1)  
        cnt = bwt_text.count(char)
        occurrences[char] = {'First Occurrence': fo, 'Count': cnt}
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
                break 
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

@app.route('/process', methods=['POST'])
def process():
    tool = request.form['tool']
    try:
        if tool == 'search':
            # Handle both text input and file upload
            if 'text_file' in request.files and request.files['text_file'].filename != '':
                file = request.files['text_file']
                if file.filename.lower().endswith(('.txt', '.fasta', '.fa')):
                    # Read file content as text
                    file_content = file.read().decode('utf-8')
                    
                    if file.filename.lower().endswith(('.fasta', '.fa')):
                        # Create a StringIO object to simulate a file object
                        fasta_file = StringIO(file_content)
                        
                        # Parse FASTA file
                        sequences = list(SeqIO.parse(fasta_file, "fasta"))
                        if not sequences:
                            return jsonify({'success': False, 'error': 'No sequences found in FASTA file'})
                        text = "".join(str(record.seq) for record in sequences)
                    else:
                        # Handle text file
                        text = file_content
                else:
                    return jsonify({'success': False, 'error': 'Invalid file type. Use .txt, .fasta, or .fa'})
            else:
                # Handle direct text input
                text = request.form['text']
            
            pattern = request.form['pattern']
            
            if not text:
                return jsonify({'success': False, 'error': 'No input text provided'})
            if not pattern:
                return jsonify({'success': False, 'error': 'No pattern provided'})
            
            occurrences = search_with_suffix_array(text, pattern)
            result = {
                'success': True,
                'occurrences': occurrences,
                'input_length': len(text),
                'matches': len(occurrences)
            }
            
            return jsonify(result)



        elif tool == 'suffix_array':
            text = request.form['text']
            sa = build_suffix_array(text)
            df = pd.DataFrame({
                'Index': sa,
                'Suffix': [text[i:] for i in sa]
            })
            return jsonify({
                'success': True,
                'table': df.to_html(classes='table table-dark', index=False)
            })

        elif tool == 'bwt':
            text = request.form['text']
            return jsonify({
                'success': True,
                'result': bwt(text)
            })

        elif tool == 'fm_index':
            text = request.form['text']
            pattern = request.form['pattern']
            fm = create_fm_index(text, pattern)
            df = pd.DataFrame([{
                'Character': k,
                'First': v['First'],
                'Count': v['Count']
            } for k, v in fm.items()])
            return jsonify({
                'success': True,
                'table': df.to_html(classes='table table-dark', index=False)
            })

        elif tool == 'dna_protein':
            dna = request.form['dna']
            protein, amino = translate_dna_to_protein(dna)
            return jsonify({
                'success': True,
                'protein': protein,
                'amino': amino
            })

        elif tool == 'alignment':
            result = perform_global_alignment(
                request.form['seq1'],
                request.form['seq2'],
                float(request.form['match']),
                float(request.form['mismatch']),
                float(request.form['open']),
                float(request.form['extend'])
            )
            return jsonify({
                'success': True,
                'alignment': result
            })

    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })

if __name__ == '__main__':
    app.run(debug=True)
