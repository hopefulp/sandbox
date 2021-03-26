#!/home/joonho/anaconda3/bin/python

from flask import Flask

app = Flask(__name__)

#@app.rout('/')
@app.route('/user/<user_name>/<int:user_id>')

#def home():
#    return 'Hello, World!'
def user(user_name, user_id):
    return f'Hello, {user_name}({user_id})!'


if __name__ == '__main__':
    app.run(debug=True)


