import unittest

from ensembl_map.util import assert_valid_position, is_ensembl_id, singleton


class TestAssertValidPosition(unittest.TestCase):
    def test_assert_valid_postion_lt_one_1(self):
        self.assertRaises(ValueError, assert_valid_position, 0)

    def test_assert_valid_postion_lt_one_2(self):
        self.assertRaises(ValueError, assert_valid_position, 0, 1)

    def test_assert_valid_postion_lt_one_3(self):
        self.assertRaises(ValueError, assert_valid_position, 1, 0)

    def test_assert_valid_postion_compare(self):
        self.assertRaises(ValueError, assert_valid_position, 100, 1)


class TestIsEnsemblId(unittest.TestCase):
    def test_is_ensembl_id_true(self):
        self.assertTrue(is_ensembl_id("ENSG00000000000"))

    def test_is_ensembl_id_true_1(self):
        self.assertTrue(is_ensembl_id("ENSG00000000000.0"))

    def test_is_ensembl_id_false(self):
        self.assertFalse(is_ensembl_id("blah"))

    def test_is_ensembl_id_false_2(self):
        self.assertFalse(is_ensembl_id("ENSG0"))


class TestSingleton(unittest.TestCase):
    @singleton
    class DummyClass:
        pass

    def test_singleton(self):
        self.assertEqual(self.DummyClass(), self.DummyClass())
