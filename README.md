Used Gemini Pro for this project. 

the project displays different proteins, but now much more substantially than in hw4. it shows tons of ligands already bound to pockets, and if users click them, API calls are made, resulting in the ligand name,
chemical formula, molecular weight, ai agent summary, and more are displayed! 
im most proud of the overall cohesive-ness of the project that feels informative to users that have an already intermediate level understanding of the pertinent biology/chemistry.
the api key (secret) is an environment variable on render.

To run locally, you should input your own api key into the app.py code (in the line w/ os.get ....)
pip install each of the requirements
change API_URL line in visualizer.html to const API_URL = "http://127.0.0.1:5000/api/molecule"; 
run app.py, open visualizer.html from local location
