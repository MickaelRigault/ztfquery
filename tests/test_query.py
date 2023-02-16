import os, logging
import unittest

from astropy.time import Time
from ztfquery import query

logging.getLogger("ztfquery.query").setLevel(logging.DEBUG)


class TestQuery(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

    def test_simple_jd_query(self):
        self.logger.info("\n\n Testing querying IRSA for time-range\n")

        zquery = query.ZTFQuery()

        jdstart = Time("2018-05-09T11:50:00").jd
        jdend = Time("2018-05-10T00:00:00").jd

        zquery.load_metadata(
            sql_query=f"seeing<2 and obsjd BETWEEN {jdstart} AND {jdend}"
        )
        n_results = len(zquery.metatable)
        n_expected = 6

        self.assertEqual(n_results, n_expected)

    def test_object_query(self):
        self.logger.info(
            "\n\n Testing querying IRSA for object coordinates within jd-ranges\n"
        )
        zquery = query.ZTFQuery()
        starttime = Time("2018-05-14").jd
        endtime = Time("2018-10-14").jd
        # Do the Query to see what exists
        zquery.load_metadata(
            radec=[276.107960, +44.130398],
            size=0.01,
            sql_query=f"fid=3 and obsjd>{starttime} and obsjd<{endtime}",
        )
        n_results = len(zquery.metatable)
        n_expected = 53

        self.assertEqual(n_results, n_expected)

    def test_ref_image_query(self):
        self.logger.info("\n\n Testing getting reference images\n")

        zquery = query.ZTFQuery()
        zquery.load_metadata(kind="ref", radec=[276.107960, +44.130398], size=0.0001)

        n_results = len(zquery.metatable)
        n_expected = 6

        self.assertEqual(n_results, n_expected)
