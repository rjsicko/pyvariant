import unittest
import os

from ensembl_map.cache import Ensembl, set_cache_dir, set_ensembl_release


# Enable after implementing test data:


# class TestSetCacheDir(unittest.TestCase):
#     def setUp(self):
#         self.old_path = os.environ.get("PYENSEMBL_CACHE_DIR")

#     def tearDown(self):
#         if self.old_path:
#             os.environ["PYENSEMBL_CACHE_DIR"] = self.old_path

#     def test_set_cache_dir(self):
#         set_cache_dir("test")
#         self.assertEqual(os.environ["PYENSEMBL_CACHE_DIR"], "test")


# class TestSetEnsemblRelease(unittest.TestCase):
#     def test_set_cache_dir(self):
#         data = set_ensembl_release(99, "human", download_if_missing=False)
#         self.assertEqual(data.release, 99)
