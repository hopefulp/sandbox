#!/home/joonho/anaconda3/bin/python

from flask import Flask

app = Flask(__name__)

@app.route('/')
#@app.route('/user/<user_name>/<int:user_id>')
#@app.route("/runVASP", method=['POST', 'GET])

#def home():
#    return 'Hello, World!'
def user(user_name, user_id):
    return f'Hello, {user_name}({user_id})!'


if __name__ == '__main__':
    #app.run(debug=True)
    app.run(host="143.248.249.73", port=5000)


