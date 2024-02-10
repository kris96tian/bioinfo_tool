# /mnt/chromeos/MyFiles/BTools
class FmIndex:
    def __init__(self, text):
        self.text = text
        self.suffix_array = self._build_suffix_array()
        self.bwt = self._construct_bwt()
        self.first_occurrence, self.count = self._create_fm_index()

    def _build_suffix_array(self):
        suffixes = [(self.text[i:], i) for i in range(len(self.text))]
        suffixes.sort()
        return [suffix[1] for suffix in suffixes]

    def _construct_bwt(self):
        return ''.join(self.text[i - 1] if i > 0 else self.text[-1]
                       for i in self.suffix_array)

    def _create_fm_index(self):
        first_occurrence = {}
        count = {}
        symbols = set(self.bwt)

        for symbol in symbols:
            first_occurrence[symbol] = self.bwt.index(symbol)
            count[symbol] = [0] * (len(self.bwt) + 1)

        for i, symbol in enumerate(self.bwt):
            for s in symbols:
                count[s][i + 1] = count[s][i] + (1 if s == symbol else 0)

        return first_occurrence, count

    def fm_substring_search(self, pattern):
        top = 0
        bottom = len(self.bwt) - 1
        for i in range(len(pattern) - 1, -1, -1):
            symbol = pattern[i]
            if symbol not in self.first_occurrence:
                return None
            top = self.first_occurrence[symbol] + self.count[symbol][top]
            bottom = self.first_occurrence[symbol] + self.count[symbol][bottom + 1] - 1
            if top > bottom:
                return None
        return (top, bottom)
