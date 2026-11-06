

LESSON1=derivace_problem.py slozeni_problem.py H2O_problem.py
LESSON2=chem_pot_1_problem.py chem_pot_2_problem.py chem_pot_3_problem.py chem_pot_4_problem.py chem_pot_5_problem.py

all: $(LESSON1) $(LESSON2)

%_problem.py: %.py
	./make_problems.py $^
