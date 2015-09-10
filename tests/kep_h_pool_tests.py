from nose.tools import *
from kep_h.kep_h_pool import *

def setup():
    pass

def teardown():
    pass

def test_empty_pool():
    pool = Pool()
    for max_length in range(5):
        assert_equal(pool.find_cycles(max_length), [])
        assert_equal(pool.find_chains(max_length), [])

def p(id):
    p.index += 1
    return Patient(id, p.index)
p.index = 0

def d(id):
    return PairedDonor(50, id)

def test_is_in():
    pool = Pool()

    p1 = p(1)
    p2 = p(2)
    d1 = d(1)
    d2 = d(2)

    pair_list = [PatientDonorPair(p1, d1),
                 PatientDonorPair(p2, d2)]
    assert p1.is_in(pair_list)
    assert p2.is_in(pair_list)
    assert d1.is_in(pair_list)
    assert d2.is_in(pair_list)

    pair_list = [PatientDonorPair(p1, d1),
                 PatientDonorPair(p1, d2)]
    assert not p2.is_in(pair_list)

    pair_list = [PatientDonorPair(p1, d1),
                 PatientDonorPair(p2, d1)]
    assert not d2.is_in(pair_list)

def test_cycle_finding():
    pool = Pool()

    p1 = p(1)
    p2 = p(2)
    p3 = p(3)
    d1 = d(1)
    d2 = d(2)
    d3 = d(3)

    for patient in [p1, p2, p3]:
        pool.patients.append(patient)

    for donor in [d1, d2, d3]:
        pool.paired_donors.append(donor)

    pool.associate_patient_with_donor(p1, d1)
    pool.associate_patient_with_donor(p2, d2)
    pool.associate_patient_with_donor(p3, d3)

    assert pool.find_cycles(3) == []
    
    for i in range(3):
        [d1, d2, d3][i-1].edges_out.append(
                DonorPatientMatch([p1, p2, p3][i], 100))
    
    assert len(pool.find_cycles(3)) == 1
    cycle = pool.find_cycles(3)[0]
    assert len(cycle.pd_pairs) == 3

    for i in range(3):
        [d1, d2, d3][i].edges_out.append(
                DonorPatientMatch([p1, p2, p3][i-1], 100))
    cycles = pool.find_cycles(3)
    assert len(cycles) == 5
    for c in cycles:
        assert len(c.pd_pairs)==c.n_transplants()
    
def test_cycle_backarcs():
    pool = Pool()

    p1 = p(1)
    p2 = p(2)
    p3 = p(3)
    d1 = d(1)
    d2 = d(2)
    d3 = d(3)

    for patient in [p1, p2, p3]:
        pool.patients.append(patient)

    for donor in [d1, d2, d3]:
        pool.paired_donors.append(donor)

    pool.associate_patient_with_donor(p1, d1)
    pool.associate_patient_with_donor(p2, d2)
    pool.associate_patient_with_donor(p3, d3)

    for i in range(3):
        [d1, d2, d3][i-1].edges_out.append(
                DonorPatientMatch([p1, p2, p3][i], 100))
    
    for i in range(3):
        [d1, d2, d3][i].edges_out.append(
                DonorPatientMatch([p1, p2, p3][i-1], 100))
    
    cycles = pool.find_cycles(3)
    for c in cycles:
        if len(c.pd_pairs) == 3:
            assert c.n_backarcs() == 3
