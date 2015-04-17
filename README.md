# boqa-negative
Adding functionality for negative phenotype annotation to Bauer et. Al's BOQA
A single diagnosis (full run through of computing each marginal given a fixed configuration of annotations) currently takes ~1.2 seconds.

To recreate results:
To generate patients:
	python generate_patient_pairs.py --out_path generated_patients_neg -N 100 --generate PATIENTS --imprecision --noise 0.5 --negative_phenotypes 
	python generate_patient_pairs.py --out_path generated_patients_neg -N 100 --generate PATIENTS --imprecision --noise 0.5

To run respective methods on each set of patients:
	submit_(k|ic_sampling|sampling|p_sampling).sh directory_with_patients output_folder

To run nofreqs on each set of patients:
	python run_net.py --patient_path directory_with_patients --out_path output_folder

hpo.py, omim.py are adapted from github.com/buske/patient-similarity
generate_patient_pairs.py is adapted from github.com/talf301/matchingsim
