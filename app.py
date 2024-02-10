from flask import Flask, request, render_template
from fm_index import FmIndex

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        text = request.form.get('text')
        substring = request.form.get('substring')

        # Input Validation (Illustrative)
        if not text:
            return render_template('index.html', error="Please enter input text.")  

        try:
            fm = FmIndex(text)
            result = fm.fm_substring_search(substring)
            if result:
                top, bottom = result
                indexes = fm.suffix_array[top:bottom + 1]
                return render_template('result.html', substring=substring, indexes=indexes, fm=fm)
            else:
                return render_template('result.html', substring=substring, indexes=None, fm=fm)

        except ValueError:  # Example of catching potential errors
            return render_template('index.html', error="An error occurred. Please check your input.")  

    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)

