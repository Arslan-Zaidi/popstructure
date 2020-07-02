
effects_i=${1}

for j in {1..20}; do \
python generate_gvalue_sib.py --chr ${j} --effect_i ${effects_i}; done
