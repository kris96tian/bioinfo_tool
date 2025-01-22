import streamlit as st
import numpy as np
from io import StringIO
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
import pandas as pd

st.set_page_config(
    page_title="Bioinformatics String Tools",
    page_icon="💾",
    layout="wide",
    initial_sidebar_state="expanded",
)

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

def main():
    st.title("Bioinformatics Tools")
    st.markdown("Generate suffix array, Burrows-Wheeler transform, and FM indices for input strings. Translate DNA sequences into their corresponding protein sequences. Perform global alignment of input strings.")

    st.sidebar.title("Select Tool")
    tool = st.sidebar.radio("Select Tool", ["Suffix Array", "Burrows-Wheeler Transform", "FM Index", "DNA to Protein", "Global Alignment", "Search with Suffix Array"])

    if tool == "Suffix Array":
        st.subheader("Suffix Array")
        text = st.text_area("Enter text", height=200)
        if st.button("Generate Suffix Array"):
            suffixes = build_suffix_array(text)
            
            data = []
            for idx in suffixes:
                suffix = text[idx:]
                data.append({
                    "Index": idx,
                    "Suffix": suffix
                })
            
            df = pd.DataFrame(data)
            st.write("Suffix Array:")
            st.table(df)

    elif tool == "Search with Suffix Array":
        st.subheader("Search with Suffix Array")
        input_option = st.radio("Input Option", ["Enter Text", "Upload File"])
        text = ""
        if input_option == "Enter Text":
            text = st.text_area("Enter text", height=200)
        else:
            uploaded_file = st.file_uploader("Upload a text file", type=["txt"])
            if uploaded_file is not None:
                text = uploaded_file.getvalue().decode("utf-8")
        
        pattern = st.text_input("Enter pattern to search")
        if st.button("Search"):
            if text and pattern:
                occurrences = search_with_suffix_array(text, pattern)
                if occurrences:
                    st.write(f"Pattern found at indices: {occurrences}")
                else:
                    st.write("Pattern not found.")
            else:
                st.write("Please provide both text and pattern.")


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
            unique_chars = []
            for char in pattern:
                if char not in unique_chars:
                    unique_chars.append(char)
            data = []
            for char in unique_chars:
                fo = occurrences.get(char, {}).get('First Occurrence', -1)
                cnt = occurrences.get(char, {}).get('Count', 0)
                data.append({
                    "Character": char,
                    "First Occurrence": fo,
                    "Count": cnt
                })
            df = pd.DataFrame(data)
            st.write("FM Index:")
            st.table(df)
            
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
        st.subheader("Needleman-Wunsch Global Alignment")
        seq1 = st.text_input("Enter the first sequence", value="AGTACG", placeholder="e.g., AGTACG")
        seq2 = st.text_input("Enter the second sequence", value="ACATAG", placeholder="e.g., ACATAG")
        
        st.write("Alignment Parameters:")
        col1, col2 = st.columns(2)
        with col1:
            match_score = st.number_input("Match Score", value=1)
            mismatch_score = st.number_input("Mismatch Score", value=-1)
        with col2:
            open_gap_score = st.number_input("Open Gap Score", value=-10)
            extend_gap_score = st.number_input("Extend Gap Score", value=-0.5)
        
        if st.button("Perform Global Alignment"):
            aligned_result = perform_global_alignment(seq1, seq2, match_score, mismatch_score, open_gap_score, extend_gap_score)
            st.code(aligned_result, language="text")


if __name__ == "__main__":
    main()
