# @kris96tian

class FmIndex:
    """
    This class represents an FM-Index, a data structure used for efficient substring
    searching in a text. It stores the Burrows-Wheeler Transform (BWT), suffix array,
    first occurrence, and count tables.

    Args:
        text (str): The input text for which the FM-Index is constructed.

    Attributes:
        text (str): The input text.
        suffix_array (list of int): The suffix array of the input text.
        bwt (str): The Burrows-Wheeler Transform of the text.
        first_occurrence (dict): A dictionary storing the first occurrence position
            of each symbol in the BWT.
        count (dict): A dictionary storing count tables for each symbol in the BWT.

    Methods:
        - build_suffix_array(): Constructs the suffix array for the given text.
        - construct_bwt(): Constructs the Burrows-Wheeler Transform (BWT) of the text.
        - create_fm_index(): Creates the first occurrence and count tables for the BWT.
        - fm_substring_search(pattern): Searches for a pattern in the BWT using FM-index.
        - display_fm_index(): Displays information about the FM-Index.

    Usage:
    fm = FmIndex("exampletext")
    fm.display_fm_index()
    result = fm.fm_substring_search("pattern")
    """

    def __init__(self, text):
        self.text = text
        self.suffix_array = self.build_suffix_array()
        self.bwt = self.construct_bwt()
        self.first_occurrence, self.count = self.create_fm_index()

    def build_suffix_array(self):
        """
        Constructs the suffix array of the input text.

        Returns:
            list of int: The suffix array.
        """
        suffixes = []
        for i in range(len(self.text)):
            suffixes.append((self.text[i:], i))
        suffixes.sort()
        suffix_array = []
        for suffix in suffixes:
            suffix_array.append(suffix[1])
        return suffix_array

    def construct_bwt(self):
        """
        Constructs the Burrows-Wheeler Transform (BWT) of the input text.

        Returns:
            str: The BWT of the text.
        """
        bwt = ''
        for i in self.suffix_array:
            if i > 0:
                bwt += self.text[i - 1]
            else:
                bwt += self.text[-1]
        return bwt

    def create_fm_index(self):
        """
        Creates the first occurrence and count tables for the BWT.

        Returns:
            dict: The first occurrence table.
            dict: The count tables for each symbol.
        """
        first_occurrence = {}
        count = {}
        symbols = set(self.bwt)
        sorted_bwt = sorted(self.bwt)
        for symbol in symbols:
            first_occurrence[symbol] = sorted_bwt.index(symbol)
        for symbol in symbols:
            count[symbol] = [0] * (len(self.bwt) + 1)
        for i, symbol in enumerate(self.bwt):  # Calc. count for each symbol
            for s in symbols:
                count[s][i + 1] = count[s][i] + 1 if s == symbol else count[s][i]
        return first_occurrence, count

    def fm_substring_search(self, pattern):
        """
        Searches for a pattern in the BWT using the FM-index.

        Args:
            pattern (str): The pattern to search for.

        Returns:
            tuple or None: A tuple (top, bottom) representing the range of occurrences
            in the BWT, or None if the pattern is not found.
        """
        top = 0
        bottom = len(self.bwt) - 1
        for i in range(len(pattern) - 1, -1, -1):
            symbol = pattern[i]
            top = self.first_occurrence[symbol] + self.count[symbol][top]
            bottom = self.first_occurrence[symbol] + self.count[symbol][bottom + 1] - 1
            if top > bottom:
                return None
        return (top, bottom)

    def display_fm_index(self):
        """
        Displays information about the FM-Index, including the suffix array, BWT,
        first occurrence, and count tables.
        """
        print("Suffix Array:", self.suffix_array)
        print("BWT:", self.bwt)
        print("First Occurrence:", self.first_occurrence)
        print("Count:")
        for symbol, counts in self.count.items():
            print(f"    {symbol}: {counts}")
