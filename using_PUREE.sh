# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# # PUREE necessitates Python <3.10 and scikit-learn <1.2.0
# Install pyenv
curl https://pyenv.run | bash
# Add to your ~/.bashrc
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(pyenv init -)"' >> ~/.bashrc
source ~/.bashrc
# Install Python 3.10
pyenv install 3.10.15

# # Preparing the venv
# Going to the right directory and venv
cd /home/eugenie-modolo/Documents/Gynet/scripts/PUREE
# Create venv with Python 3.10
pyenv local 3.10.15 
# to check the python version: 
# python --version
python -m venv puree_env
# to deactivate the venv: 
# deactivate
# to delete the venv: 
# rm -rf puree_env
source puree_env/bin/activate
pip install -r requirements.txt

# # Running PUREE
# Going to the right directory and venv
cd /home/eugenie-modolo/Documents/Gynet/scripts/PUREE
# activate the venv
source puree_env/bin/activate
# run from inside the installation directory
python3 predict_purity.py --data_path ./counts_Gynet_puree.csv \
			  			  --output puree_result_Gynet.tsv 
python3 predict_purity.py --data_path ./counts_Lapnet_puree.csv \
						  --output puree_result_Lapnet.tsv
python3 predict_purity.py --data_path ./counts_control_puree.csv \
						  --output puree_result_control.tsv
python3 predict_purity.py --data_path /home/eugenie-modolo/Documents/Donnees_Netris/scripts/counts_IHCAF_puree.csv \ 
                          --output /home/eugenie-modolo/Documents/Donnees_Netris/scripts/puree_result_IHCAF.tsv
python3 predict_purity.py --data_path /home/eugenie-modolo/Documents/Donnees_Netris/scripts/counts_BCOR_puree.csv \ 
                          --output /home/eugenie-modolo/Documents/Donnees_Netris/scripts/puree_result_BCOR.tsv
# then deactivate the venv
deactivate
