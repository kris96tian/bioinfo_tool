import streamlit as st
import numpy as np
from io import StringIO
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align

st.set_page_config(
    page_title="Bioinformatics String Tools",
    page_icon="💾",
    layout="wide",
    initial_sidebar_state="expanded",
)

def perform_global_alignment(seq1, seq2):
    """
    Performs a global alignment between the two input sequences.
    Returns the aligned sequences as a string.
    """
    # Perform the global alignment
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(str(seq1), str(seq2))
    alignment = alignments[0]
    print(alignment)
    return str(alignment)

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
    
    occurrences = {}
    
    for char in pattern:
        first_occurrence = -1
        count = 0
        left, right = 0, len(text) - 1
        
        while left <= right:
            mid = (left + right) // 2
            if len(text[suffix_array[mid]:]) >= len(char):
                if text[suffix_array[mid]:suffix_array[mid] + len(char)] == char:
                    # Find the leftmost occurrence
                    left_idx = mid
                    while left_idx > 0 and text[suffix_array[left_idx-1]:suffix_array[left_idx-1] + len(char)] == char:
                        left_idx -= 1
                    first_occurrence = suffix_array[left_idx]
                    
                    # Count the number of occurrences
                    right_idx = left_idx
                    while right_idx < len(text) and text[suffix_array[right_idx]:suffix_array[right_idx] + len(char)] == char:
                        right_idx += 1
                    count = right_idx - left_idx
                    break
                elif text[suffix_array[mid]:suffix_array[mid] + len(char)] < char:
                    left = mid + 1
                else:
                    right = mid - 1
            elif text[suffix_array[mid]:] < char:
                left = mid + 1
            else:
                right = mid - 1
        
        occurrences[char] = (first_occurrence, count)
    
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

def main():
    st.title("Bioinformatics Tools")
    st.markdown("Generate suffix array, Burrows-Wheeler transform, and FM indices for input strings. Translate DNA sequences into their corresponding protein sequences. Perform global alignment of input strings.")

    st.sidebar.title("Select Tool")
    tool = st.sidebar.radio("Select Tool", ["Suffix Array", "Burrows-Wheeler Transform", "FM Index", "DNA to Protein", "Global Alignment"])

    if tool == "Suffix Array":
        st.subheader("Suffix Array")
        text = st.text_area("Enter text", height=200)
        if st.button("Generate Suffix Array"):
            suffixes = build_suffix_array(text)
            st.write("Suffix Array:")
            st.write(suffixes)

    elif tool == "Burrows-Wheeler Transform":
        st.subheader("Burrows-Wheeler Transform")
        text = st.text_area("Enter text", height=200)
        if st.button("Generate BWT"):
            bwt_text = bwt(text)
            st.write("Burrows-Wheeler Transform:")
            st.write(bwt_text)

    elif tool == "FM Index":
        st.subheader("FM Index")
        text = st.text_area("Enter text", height=200)
        pattern = st.text_input("Enter pattern")
        if st.button("Generate FM Index"):
            occurrences = create_fm_index(text, pattern)
            result = "FM Index:\n"
            result += "{\n"
            for char in pattern:
                if char in occurrences:
                    first_occurrence, count = occurrences[char]
                    result += f" \n   '{char}': {{'First Occurrence': {first_occurrence}, 'Count': {count}}},\n"
                else:
                    result += f" \n   '{char}': {{'First Occurrence': -1, 'Count': 0}},\n"
            result += "}\n"
            st.write(result)

    elif tool == "DNA to Protein":
        st.subheader("DNA to Protein Translation")
        dna_sequence = st.text_area("Enter DNA sequence", height=200)
        if st.button("Translate to Protein"):
            protein_sequence, amino_acid_sequence = translate_dna_to_protein(dna_sequence)
            st.write("Protein Sequence:")
            st.write(protein_sequence)
            st.write("Amino Acid Sequence:")
            st.write(amino_acid_sequence)

    elif tool == "Global Alignment":
        st.subheader("Needleman-Wunsch Alignment")
        seq1 = st.text_input("Enter the first sequence")
        seq2 = st.text_input("Enter the second sequence")
        if st.button("Perform Global Alignment"):
            aligned_sequence = perform_global_alignment(seq1, seq2)
            st.write(aligned_sequence)

if __name__ == "__main__":
    main()
