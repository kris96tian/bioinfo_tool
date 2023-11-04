from flask import Flask, request, render_template
from fm_index import FmIndex

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        text = request.form.get('text')
        substring = request.form.get('substring')
        fm = FmIndex(text)
        result = fm.fm_substring_search(substring)
        if result:
            top, bottom = result
            indexes = fm.suffix_array[top:bottom + 1]
            return render_template('result.html', substring=substring, indexes=indexes, fm=fm)
        else:
            return render_template('result.html', substring=substring, indexes=None, fm=fm)
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
