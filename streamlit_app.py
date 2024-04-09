import streamlit as st
import subprocess
import os


st.set_page_config(
    page_title="Bioinformatics String Tools",
    page_icon="💾",
    layout="wide",
    initial_sidebar_state="expanded",
)

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
    st.markdown("Generate suffix array, Burrows-Wheeler transform, and FM indices for input strings. You can also translate DNA sequences into their corresponding proteine and aminoaci sequences, as well as get the global alignments with scores of two fasta files.")

    # Sidebar with sections for each tool
    tool = st.sidebar.radio("Select Tool", ["Suffix Array", "Burrows-Wheeler Transform", "FM Index", "DNA to Protein", "Sequence Alignment"])

    # Suffix Array Tool
    if tool == "Suffix Array":
        st.header("Get suffix array")
        text_suffix = st.text_area("Enter text", height=200)
        if st.button("Generate Suffix Array"):
            suffixes = build_suffix_array(text_suffix)
            st.write("Suffix Array:")
            st.write(suffixes)

    # Burrows-Wheeler Transform Tool
    elif tool == "Burrows-Wheeler Transform":
        st.header("Get the BWT of string:")
        text = st.text_area("Enter text", height=200)
        if st.button("Generate BWT"):
            bwt_text = bwt(text)
            st.write("Burrows-Wheeler Transform:")
            st.write(bwt_text)

    # FM Index Tool
    elif tool == "FM Index":
        st.header("Generate the FM Index:")
        text_fm = st.text_area("Enter text", height=200)
        pattern = st.text_input("Enter pattern")
        if st.button("0 : First Occurrence, 1 : Count"):
            occurrences = create_fm_index(text_fm, pattern)
            st.write("FM Index:")
            st.write(occurrences)

    # DNA to Protein Tool
    elif tool == "DNA to Protein":
        st.header("Translate DNA Sequence:")
        dna_sequence = st.text_area("Enter DNA sequence", height=200)
        if st.button("Translate to Protein"):
            protein_sequence, amino_acid_sequence = translate_dna_to_protein(dna_sequence)
            st.write("Protein Sequence:")
            st.write(protein_sequence)
            st.write("Amino Acid Sequence:")
            st.write(amino_acid_sequence)

    # Sequence Alignment Tool
    elif tool == "Sequence Alignment":
        st.title("Get Sequence Alignment:")
        file1 = st.file_uploader("Upload File 1", type=["txt"])
        file2 = st.file_uploader("Upload File 2", type=["txt"])
        match = st.number_input("Match Score", value=3, step=1)
        mismatch = st.number_input("Mismatch Score", value=-1, step=1)
        gap = st.number_input("Gap Score", value=-2, step=1)
        if st.button("Align"):
            if file1 is not None and file2 is not None:
                file1.seek(0)
                file2.seek(0)
                with open("file1.txt", "wb") as f:
                    f.write(file1.read())
                with open("file2.txt", "wb") as f:
                    f.write(file2.read())
                try:
                    result = subprocess.run(['./align', 'file1.txt', 'file2.txt', str(int(match)), str(int(mismatch)), str(int(gap))], capture_output=True, text=True)
                    os.remove('file1.txt')
                    os.remove('file2.txt')
                    st.write("Click the button below to download the alignment results.")
                    st.download_button(
                        "Download Results",
                        result.stdout,
                        "alignment_results.txt",
                        "text/plain",
                    )
                except subprocess.CalledProcessError as e:
                    st.error(f"Error: {e.stderr}")
            else:
                st.warning("Please upload both input files.")

    # Footer
    st.write("\n\n\n")
    st.write("Created by Kristian Alikaj")

if __name__ == '__main__':
    main()
