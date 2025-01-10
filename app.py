from flask import Flask, render_template, request, redirect, url_for, Response

# import functions from science api: search and show images
from similarity import findSimilarCompounds, findSimilarity, smiles_to_svg

app = Flask(__name__)

# show main search page
@app.route("/", methods=["GET"])
def root():
    # user input default to none
    user_input = None

    # check if user input is present
    # unsure if it's a POST or GET
    if request.method == "POST" and '' in request.form:
        # display name
        smiles_str = request.form['user_str']    # change input name if updated
        # if so, show all matching compounds
        if smiles_str:
            try:
                # get the matching ones
                displayed_comp = 

                # update for displaying on page
                user_input = smiles_str

            except ValueError:
                user_input = "Invalid input!"

        return render_template("result.html", user_input=user_input, compounds=displayed_comp)     # include list of compounds


    # if no user input, display home page
    return render_template("base.html")


# render page to show results
@app.route("/similar-comp", methods=['GET'])
def show_result(mol_name):
    #
    return render_template('.html', )   # update html name and include compound info in a dict

# for getting comp image
@app.route("/comp-image", methods=["GET"])
def get_image():
    return Response(smiles_to_svg(), mimetype="image/svg+xml")