# kep-hierarchical

This program aims to re-implement the hierarchical model for kidney-exchange
optimisation described in _Paired and Altruistic
Kidney Donation in the UK: Algorithms and Experimentation_. David Manlove
and Gregg O'Malley. J. Exp. Algorithmics 19, Article 2.6 (January 2015)
[Link](http://dl.acm.org/citation.cfm?doid=2627368.2670129).

*Please note that the program in this repository is at an early stage of
development and is intented for experimental work only.*

## Objective criteria

Currently, the five objective criteria used in the UK National Living
Donor Kidney Sharing Scheme are implemented. In the descriptions below,
"3-cycles" includes chains of length 2.

- *effective*: maximise the number of effective two-cycles
- *size*: maximise the number of transplants (including transplants
  from unused altruistic donors)
- *3way*: minimise the number of 3-cycles
- *backarc*: maximise the number of backarcs in 3-cycles
- *weight*: maximise the overall weight of the solution, using the
  NLDKSS weight function

It should be reasonably straightforward to add new criteria to
`kep_h/optimality_criteria.py`, if required, by adding a name to
the get\_criterion function and adding a new subclass of OptCriterion.

## Command-line usage

The program currently only accepts JSON input files in the format
described [here](http://kidney.optimalmatching.com/api/input_format).

For help, type:
```
python kep_h.py -h
```

An example invocation, with maximum cycle length 3, maximum chain
length 1 (so-called "short chains"), stopping after 10 optimal solutions
have been found:
```
python kep_h.py -f input-file.json -c effective:size:3way -e 3 -n 1 -m 10
```

The output will look something like this:
```
[1(100) 2(200)]
[(1000) 3(300) 4(400)]
altruist 2000 
```

In the first line of this example, the patient-donor pair with patient 1
and donor 100 is selected in a 2-cycle with the pair (2, 200). In the second
line, the altruistic donor 1000 initiates a 2-chain with the patient-donor
pairs (3, 300) and (4, 400). The final line indicates that the altruistic
donor 2000 is not selected for a chain.
