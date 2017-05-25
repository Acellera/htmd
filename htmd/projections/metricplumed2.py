# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.projections.projection import Projection
from htmd.molecule.molecule import Molecule
from htmd.util import ensurelist

from abc import ABC
import logging
import numpy
import subprocess
import shutil
import os
import tempfile

logger = logging.getLogger(__name__)

cvcounter = 1


# Local utility functions --------------------------------------------------------

def _getTempDirName(prefix=""):
    return os.path.join(tempfile._get_default_tempdir(),
                        prefix + next(tempfile._get_candidate_names()))


def _getPlumedRoot():
    """ Return path to plumed executable, or raise an exception if not found. """
    return _getPlumedInfo("root")


def _getPlumedInfo(what):
    """ Return selected info from plumed executable, or raise an exception if not found. """
    info = subprocess.check_output(["plumed", "--standalone-executable", "info", "--" + what])
    info = info.strip().decode("utf-8")
    return info

def _printDFS(n):
    """ Depth-first traversal of nodes """
    # https://en.wikipedia.org/wiki/Topological_sorting#Depth-first_search
    _L = []
    _tempmarked = dict()
    _marked = dict()

    def _printDFS_aux(nn):
        nonlocal _tempmarked,_marked,_L
        if nn in _tempmarked:
            raise Exception("Not a DAG! There may be cyclic dependencies in the colvars.")
        if not nn in _marked:
            _tempmarked[nn]=True
            if hasattr(nn,'prereq'):
                for d in nn.prereq:
                    _printDFS_aux(d)
            _marked[nn]=True
            del _tempmarked[nn]
            _L = [str(nn)] + _L

    _printDFS_aux(n)
    return list(reversed(_L))


# Static utility functions ---------------------------------------------------------
def genTemplate(action, include_optional=False):
        """ Return the template for the given action

        Parameters
        ----------
        action : str
            The action to be documented
        include_optional : bool
            Whether to include optional arguments

        Examples
        --------
        >>> genTemplate("GYRATION")
        'GYRATION ATOMS=<atom selection> TYPE=RADIUS'
        """

        cl = ["plumed", "--standalone-executable", "gentemplate", "--action", action]
        if include_optional:
            cl.append("--include-optional")
        info = subprocess.check_output(cl)
        info = info.strip().decode("utf-8")
        return info


def manual(action):
    """ Return the manual for the given action.
    
    The manual is returned as a pseudo-HTML string. 
    
    Bugs: prints some text on stderr. The returned string should be reformatted to a more readable format.

    Parameters
    ----------
    action : str
        The action to be documented

    Examples
    --------
    >> manual("GYRATION")
    """

    cl = ["plumed", "--standalone-executable", "manual", "--action", action]
    info = subprocess.check_output(cl)
    info = info.strip().decode("utf-8")
    return info


# Plumed statement wrappers --------------------------------------------------------

# ABC should prevent instantiation but doesn't
class PlumedStatement(ABC):
    """ Abstract base class for Plumed statements. Do not use directly. """
    def __init__(self):
        self.prereq = []

    def __str__(self):
        return "# Rendered PlumedStatement"


class PlumedCV(PlumedStatement):
    """ Define a Plumed2 CV.

    The arguments, optional and mandatory, are passed as python named parameters. The argument
    values can be of type (see examples):

     - string or int   (passed as they are)
     - bool            (the keyword is enabled, with no value (eg: PBC=True becomes PBC))
     - PlumedGroups    (passed by label, prepended to the CV definition)
     - PlumedCOM       (passed by label, prepended to the CV definition)
     - Molecule        (converted into a PlumedGroup)
     - list containing any of the above (each element is converted as above, then they are listed separated by comma)

    Parameters
    ----------
    cv : str
        The CV action, as a string (e.g.: "DISTANCE"). (PLUMED is Case-insensitive.)
    label : str
        The label assigned to the CV
    args :
        Named arguments and keywords to the CV (see details).
    verbatim : str, optional
        Code which will be added as-is to the CV line

    Examples
    --------
    >>> PlumedCV("GYRATION", "rgyr", ATOMS="10-20", TYPE="RADIUS")   # As string
    rgyr: GYRATION ATOMS=10-20 TYPE=RADIUS
    >>> m=Molecule("3ptb")
    >>> grp=PlumedGroup(m,"grp","serial 10 to 20")
    >>> PlumedCV("GYRATION", "rgyr2", ATOMS=grp, TYPE="ASPHERICITY", NOPBC=True)  # As PlumedGroup
    grp: GROUP ATOMS=10,11,12,13,14,15,16,17,18,19,20
    rgyr2: GYRATION ATOMS=grp NOPBC TYPE=ASPHERICITY
    >>> protCA=PlumedCOM(m,"protCA","chain A and name CA")
    >>> lig=PlumedCOM(m,"lig","resname BEN and noh")
    >>> PlumedCV("DISTANCE", "dist", ATOMS=[protCA,lig])  # List of groups
    lig: COM ATOMS=1632,1633,1634,1635,1636,1637,1638,1639,1640
    protCA: COM ATOMS=2,10,17,21,25,37,44,50,54,59,67,74,81,88,100,109,116,122,130,138,144,148,160,170,181,187,191,195,201,209,217,225,231,240,254,261,268,274,279,284,294,300,312,321,327,331,339,348,355,366,374,378,387,395,403,411,419,426,433,442,446,454,463,472,483,491,497,502,508,517,523,531,538,548,555,561,573,581,587,595,602,610,618,626,634,642,650,658,666,675,683,692,698,703,708,714,722,730,736,747,754,759,765,773,779,787,794,801,807,813,818,824,829,833,840,849,855,863,871,877,881,895,899,907,914,923,929,935,939,946,952,964,971,979,986,994,1003,1009,1017,1026,1031,1038,1046,1054,1060,1068,1074,1080,1086,1095,1101,1106,1118,1125,1129,1138,1146,1153,1159,1167,1175,1186,1192,1197,1201,1213,1221,1230,1234,1238,1247,1255,1261,1267,1276,1280,1288,1294,1298,1302,1309,1316,1323,1329,1335,1339,1348,1356,1365,1369,1377,1384,1390,1404,1408,1414,1418,1424,1429,1438,1447,1455,1464,1471,1475,1482,1494,1501,1510,1517,1523,1531,1543,1550,1556,1570,1578,1587,1596,1603,1611,1616,1622,1631
    dist: DISTANCE ATOMS=protCA,lig
    >>> PlumedCV("GYRATION", "rgyr3", ATOMS=m)            # Convert Molecule implicitly
    lab_1: GROUP ATOMS=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,730,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778,779,780,781,782,783,784,785,786,787,788,789,790,791,792,793,794,795,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,829,830,831,832,833,834,835,836,837,838,839,840,841,842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,857,858,859,860,861,862,863,864,865,866,867,868,869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890,891,892,893,894,895,896,897,898,899,900,901,902,903,904,905,906,907,908,909,910,911,912,913,914,915,916,917,918,919,920,921,922,923,924,925,926,927,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,959,960,961,962,963,964,965,966,967,968,969,970,971,972,973,974,975,976,977,978,979,980,981,982,983,984,985,986,987,988,989,990,991,992,993,994,995,996,997,998,999,1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,1043,1044,1045,1046,1047,1048,1049,1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1060,1061,1062,1063,1064,1065,1066,1067,1068,1069,1070,1071,1072,1073,1074,1075,1076,1077,1078,1079,1080,1081,1082,1083,1084,1085,1086,1087,1088,1089,1090,1091,1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1111,1112,1113,1114,1115,1116,1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1127,1128,1129,1130,1131,1132,1133,1134,1135,1136,1137,1138,1139,1140,1141,1142,1143,1144,1145,1146,1147,1148,1149,1150,1151,1152,1153,1154,1155,1156,1157,1158,1159,1160,1161,1162,1163,1164,1165,1166,1167,1168,1169,1170,1171,1172,1173,1174,1175,1176,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1213,1214,1215,1216,1217,1218,1219,1220,1221,1222,1223,1224,1225,1226,1227,1228,1229,1230,1231,1232,1233,1234,1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1272,1273,1274,1275,1276,1277,1278,1279,1280,1281,1282,1283,1284,1285,1286,1287,1288,1289,1290,1291,1292,1293,1294,1295,1296,1297,1298,1299,1300,1301,1302,1303,1304,1305,1306,1307,1308,1309,1310,1311,1312,1313,1314,1315,1316,1317,1318,1319,1320,1321,1322,1323,1324,1325,1326,1327,1328,1329,1330,1331,1332,1333,1334,1335,1336,1337,1338,1339,1340,1341,1342,1343,1344,1345,1346,1347,1348,1349,1350,1351,1352,1353,1354,1355,1356,1357,1358,1359,1360,1361,1362,1363,1364,1365,1366,1367,1368,1369,1370,1371,1372,1373,1374,1375,1376,1377,1378,1379,1380,1381,1382,1383,1384,1385,1386,1387,1388,1389,1390,1391,1392,1393,1394,1395,1396,1397,1398,1399,1400,1401,1402,1403,1404,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1432,1433,1434,1435,1436,1437,1438,1439,1440,1441,1442,1443,1444,1445,1446,1447,1448,1449,1450,1451,1452,1453,1454,1455,1456,1457,1458,1459,1460,1461,1462,1463,1464,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1492,1493,1494,1495,1496,1497,1498,1499,1500,1501,1502,1503,1504,1505,1506,1507,1508,1509,1510,1511,1512,1513,1514,1515,1516,1517,1518,1519,1520,1521,1522,1523,1524,1525,1526,1527,1528,1529,1530,1531,1532,1533,1534,1535,1536,1537,1538,1539,1540,1541,1542,1543,1544,1545,1546,1547,1548,1549,1550,1551,1552,1553,1554,1555,1556,1557,1558,1559,1560,1561,1562,1563,1564,1565,1566,1567,1568,1569,1570,1571,1572,1573,1574,1575,1576,1577,1578,1579,1580,1581,1582,1583,1584,1585,1586,1587,1588,1589,1590,1591,1592,1593,1594,1595,1596,1597,1598,1599,1600,1601,1602,1603,1604,1605,1606,1607,1608,1609,1610,1611,1612,1613,1614,1615,1616,1617,1618,1619,1620,1621,1622,1623,1624,1625,1626,1627,1628,1629,1631,1632,1633,1634,1635,1636,1637,1638,1639,1640,1641,1642,1643,1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,1660,1661,1662,1663,1664,1665,1666,1667,1668,1669,1670,1671,1672,1673,1674,1675,1676,1677,1678,1679,1680,1681,1682,1683,1684,1685,1686,1687,1688,1689,1690,1691,1692,1693,1694,1695,1696,1697,1698,1699,1700,1701,1702
    rgyr3: GYRATION ATOMS=lab_1
    """

    def __init__(self, cv, label, verbatim=None, **kw):
        self.label = label
        self.cv = cv
        self.args = kw
        self.verbatim = verbatim
        self.prereq = []

        for k in self.args:
            v = self.args[k]
            # Check all possible RHS types...
            # If the arg is a Plumed group, add it to prereq and replace by label
            if isinstance(v, PlumedGenericGroup):
                self.prereq.append(v)
                self.args[k] = v.label
            # If the args is a Molecule object, convert to a PlumedGroup
            elif isinstance(v, Molecule):
                tmpGrp = PlumedGroup(label=None, mol=v, sel="all")
                self.prereq.append(tmpGrp)
                self.args[k] = tmpGrp.label
            # Boolean, as in PBC=True, is ok
            elif isinstance(v, bool):
                pass
            # String is ok
            elif isinstance(v, str):
                pass
            elif isinstance(v, int):
                self.args[k] = str(v)
            # Ditto if it is a list-like object, plus expand with commas
            elif hasattr(v, '__iter__'):
                for l in range(len(v)):
                    le = v[l]
                    if isinstance(le, PlumedGenericGroup):
                        self.prereq.append(le)
                        v[l] = le.label
                    elif isinstance(le, PlumedCV):
                        # e.g. for COMBINE
                        self.prereq.append(le)
                        if hasattr(le,"prereq"):
                            self.prereq = le.prereq+self.prereq
                        v[l] = le.label
                    elif isinstance(le, Molecule):
                        tmpGrp = PlumedGroup(label=None, mol=le, sel="all")
                        self.prereq.append(tmpGrp)
                        v[l] = tmpGrp.label
                    elif isinstance(le, str):
                        pass  # already a string
                    elif isinstance(le, int):
                        v[l] = str(le)
                    else:
                        raise TypeError("Unexpected type {} passed to argument {} (list context)".
                                        format(str(type(le)),k))
                self.args[k] = ",".join(v)
            else:
                raise TypeError("Unexpected type {} passed to argument {}".
                                        format(str(type(v)),k))

    def __str__(self):
        r = ""

        # Label
        r = r + self.label + ": " + self.cv + " "

        # Code
        if self.verbatim:
            r = r + self.verbatim + " "

        # Args
        for k, v in sorted(self.args.items()):
            if v == True:
                r = r + k + " "
            else:
                r = r + k + "=" + str(v) + " "

        return r.strip()


    def genTemplate(self, include_optional=False):
        """ Return the template for the given action

        Examples
        --------
        >>> rg=PlumedCV("GYRATION", "rgyr3", ATOMS=[1,2,3])
        >>> rg.genTemplate()
        'GYRATION ATOMS=<atom selection> TYPE=RADIUS'
        """

        tmpl = genTemplate(self.cv,include_optional=include_optional)
        return tmpl


class PlumedMolinfo(PlumedStatement):
    """ Add a MOLINFO statement. 
    
    The supplied Molecule will also be written to a temporary file, and deleted 
    when the object goes off scope.
    
    Note: this may not behave well under multiprocessing. 
    
    Parameters
    ----------
    mol: Molecule
        Will be written to file, and the corresponding MOLINFO statement generated.
        
    Examples
    --------
    >>> m = Molecule("1kdx")
    >>> mc = m.copy()
    >>> mc.dropFrames(keep=[0])
    >>> mc.set("resname","SER","resname SEP")
    >>> molinfo = PlumedMolinfo(mc)
    >>> ah = PlumedCV("ALPHARMSD",RESIDUES="119-146", R_0="1.0", label="ah")
    >>> ah_metric = MetricPlumed2([molinfo,ah])
    >>> print(ah_metric)
    MOLINFO STRUCTURE=/var/folders/qz/7p0f8wdj4zdd8nwxm89xzhy80000gn/T/tmpqqhsi644.pdb
    ah: ALPHARMSD RESIDUES=all R_0=1.0
    # Rendered PlumedStatement
    >>> ah_metric.project(m)        # doctest: +ELLIPSIS
    array([[ 20.10917664],
           [ 19.43626404],
           [ 20.09723854], ...
    """

    def __init__(self, mol):
        self.localmol = mol.copy()

    def __str__(self):
        # TODO Need to find a better pickle-compatible way which is deleted at end
        pdbfp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        self.localmol.write(pdbfp.name)
        return "MOLINFO STRUCTURE={}".format(pdbfp.name)


class PlumedVerbatim(PlumedStatement):
    """ An arbitrary Plumed statement, as a string. 
    
    Parameters
    ----------
    txt:
        The statement
    comment: bool
        Whether it is a comment.
    """
    def __init__(self, txt, comment=False):
        if comment:
            self.txt = "# "+txt
        else:
            self.txt = txt

    def __str__(self):
        return self.txt


class PlumedGenericGroup(PlumedStatement):
    """ Abstract class from which PLUMED groups are inherited. Do not use directly. """

    def __init__(self, mol, label, sel, type=""):
        global cvcounter
        al = mol.get("serial", sel)
        al = list(al)
        if label:
            self.label = label
        else:
            self.label = "lab_" + str(cvcounter)
            cvcounter = cvcounter + 1
        self.mol = mol
        self.sel = sel
        self.code = "%s: %s ATOMS=%s" % (self.label, type, ",".join(map(str, al)))

    def __str__(self):
        return self.code


class PlumedGroup(PlumedGenericGroup):
    """ An atom GROUP for use in the Plumed interface

    Parameters
    ----------
    mol: Molecule
        The molecule
    label: str
        The label assigned to the group. Autogenerated if None.
    sel: str, optional
        The atom selection defining the group (defaults to the whole molecule)

    Example
    -------
    >>> m=Molecule("3PTB")
    >>> PlumedGroup(m,"ben","resname BEN")
    ben: GROUP ATOMS=1632,1633,1634,1635,1636,1637,1638,1639,1640
    >>> g=PlumedGroup(m, label=None)   # Autogenerate label
    >>> g.label
    'lab_2'
    """

    def __init__(self, mol, label, sel="all"):
        return super(PlumedGroup, self).__init__(mol, label, sel, "GROUP")


class PlumedCOM(PlumedGenericGroup):
    """ An atom center-of-mass for use in the Plumed interface

    Parameters
    ----------
    mol: Molecule
        The molecule
    label: str
        The label assigned to the group
    sel: str, optional
        The atom selection defining the group

    Example
    -------
    >>> m=Molecule("3PTB")
    >>> PlumedCOM(m,"ben_cm","resname BEN")
    ben_cm: COM ATOMS=1632,1633,1634,1635,1636,1637,1638,1639,1640
    """

    def __init__(self, mol, label, sel="all"):
        return super(PlumedCOM, self).__init__(mol, label, sel, "COM")


# Plumed projector --------------------------------------------------------

class MetricPlumed2(Projection):
    """ Calculates generic collective variables through Plumed 2

    The collective variables are defined in PLUMED 2's syntax. PLUMED needs be installed
    separately; see http://www.plumed.org/.

    The script can be defined as
     * a string, or a list of strings (which are concatenated), or
     * a list of PlumedCV objects (see PlumedCV for examples)
     
     TODO: allow selection of CVs to print

    Parameters
    ----------
    plumed_inp :
        The PLUMED script defining CVs - string, list of strings or list of PlumedCV objects

    Examples
    --------
    >>> dd = htmd.home(dataDir="adaptive")
    >>> fsims = htmd.simlist([dd + '/data/e1s1_1/', dd + '/data/e1s2_1/'], dd + '/generators/1/structure.pdb')
    >>> metr = Metric(fsims)
    >>> metr.set(MetricPlumed2( ['d1: DISTANCE ATOMS=2,3', 'd2: DISTANCE ATOMS=5,6'])) # As strings
    >>> data=metr.project()
    >>> data.dat
    array([ array([[ 1.68597198,  1.09485197], ...
    """

    def __init__(self, plumed_inp):
        # I am not sure at all about opening files here is good style
        self._precalculation_enabled = False
        self._plumed_exe = shutil.which("plumed")
        self.colvar = None
        self.cvnames = None
        self.stmt = None

        try:
            pp = _getPlumedRoot()
            logger.info("Plumed path is " + pp)
        except Exception as e:
            raise Exception("To use MetricPlumed2 please ensure PLUMED 2's executable is installed and in path")

        # Sanitize if single element
        if type(plumed_inp) == str:
            self._plumed_inp = plumed_inp
        else:
            # This should keep the CVs etc in scope
            self.stmt = PlumedStatement()
            self.stmt.prereq = ensurelist(plumed_inp)
            stmts = _printDFS(self.stmt)
            self._plumed_inp = "\n".join(stmts)

    def __str__(self):
        return self._plumed_inp

    def _readColvar(self):
        # Assumptions: file begins with #! FIELDS time
        # first line is time
        # only colvars follow
        assert self.colvar, "colvar variable not defined"
        data = []
        with open(self.colvar, "r") as file:
            headerline = file.readline()
            self.cvnames = headerline.strip().replace('#! FIELDS time ', '').split()

            for line in file:
                if line[0] == "#":
                    continue
                cols_str = line.split()[1:]
                cols = [float(x) for x in cols_str]
                data.append(cols)
        return numpy.array(data, dtype=numpy.float32)

    # Only called if single topology
    def _precalculate(self, mol):
        logger.info("In _precalculate")
        self._precalculation_enabled = True

    def getMapping(self, mol):
        """ Return the labels of the colvars used in this projection.

        Can only be used after the projection has been executed.
        
        TODO: Update to pandas.

        Returns
        -------
        cvnames
            A list of cv names
        """
        if self.cvnames:
            return self.cvnames
        else:
            logger.warning("MetricPlumed's getMapping can only be called after the projection")
            # raise Exception("MetricPlumed's getMapping can only be called after the projection")

    # Arguments are actually self, mol
    def project(self, mol, debug=False):
        """ Project molecule.

        Parameters
        ------------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>`
            A :class:`Molecule <htmd.molecule.molecule.Molecule>` object to project.
        debug : bool
            Do not delete intermediate files.

        Returns
        -------
        data : np.ndarray
            An array containing the projected data.
        """

        logger.debug("_precalculate was called? %d" % self._precalculation_enabled)

        # --standalone-executable driver --box 100000,100000,100000 --mf_dcd /var/tmp/vmdplumed.8003/temp.dcd
        # --pdb /var/tmp/vmdplumed.8003/temp.pdb --plumed /var/tmp/vmdplumed.8003/META_INP

        td = _getTempDirName("metricplumed2-")
        os.mkdir(td)

        # PDB
        pdb = os.path.join(td, "temp.pdb")
        mol.write(pdb)

        # XTC
        xtc = os.path.join(td, "temp.xtc")
        mol.write(xtc)
        logger.debug("Done writing %d frames in %s" % (mol.numFrames, xtc))

        # Colvar
        colvar = os.path.join(td, "temp.colvar")
        self.colvar = colvar
        logger.debug("Colvar file is " + colvar)

        # Metainp
        metainp = os.path.join(td, "temp.metainp")
        metainp_fp = open(metainp, "w+")
        metainp_fp.write("UNITS  LENGTH=A  ENERGY=kcal/mol  TIME=ps\n")
        metainp_fp.write(self._plumed_inp)
        metainp_fp.write('\nPRINT ARG=* FILE=%s\n# FLUSH STRIDE=1\n' % colvar)
        metainp_fp.close()

        cmd = [self._plumed_exe, '--standalone-executable',
               'driver',
               '--mf_xtc', xtc,
               '--pdb', pdb,
               '--plumed', metainp]
        logger.debug("Invoking " + " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            logger.error("Error from PLUMED (stdout): " + e.stdout.decode("utf-8"))
            logger.error("Error from PLUMED (stderr):" + e.stderr.decode("utf-8"))
            logger.error("Leaving temporary data in " + td)
            errstr=[x for x in e.stdout.decode("utf-8").split("\n") if "PLUMED: ERROR" in x]
            raise Exception("".join(errstr))

        data = self._readColvar()
        if debug:
            logger.info("Leaving temporary data in " + td)
        else:
            shutil.rmtree(td)

        return data


# Test --------------------------------------------------------

if __name__ == "__main__":
    import sys
    import numpy as np
    from htmd import *
    from htmd.projections.metricplumed2 import *

    try:
        _getPlumedRoot()
    except:
        print("Tests in %s skipped because plumed executable not found." % __file__)
        sys.exit()



    # Simlist
    dd = htmd.home(dataDir="adaptive")
    fsims = htmd.simlist([dd + '/data/e1s1_1/', dd + '/data/e1s2_1/'],
                         dd + '/generators/1/structure.pdb')
    metr = Metric(fsims)
    metr.set(MetricPlumed2(
        ['d1: DISTANCE ATOMS=2,3',
         'd2: DISTANCE ATOMS=5,6']))
    data2 = metr.project()


    # One simulation
    testpath=os.path.join(htmd.home(), 'data', '1kdx')
    mol = Molecule(os.path.join(testpath, '1kdx_0.pdb'))
    mol.read(os.path.join(htmd.home(), 'data', '1kdx', '1kdx.dcd'))

    metric = MetricPlumed2(['d1: DISTANCE ATOMS=1,200',
                            'd2: DISTANCE ATOMS=5,6'])
    data = metric.project(mol)
    ref = np.array([0.536674, 21.722393, 22.689391, 18.402114, 23.431387, 23.13392, 19.16376, 20.393544,
                    23.665517, 22.298349, 22.659769, 22.667669, 22.484084, 20.893447, 18.791701,
                    21.833056, 19.901318])
    assert np.all(np.abs(ref - data[:, 0]) < 0.01), 'Plumed demo calculation is broken'


    #    metric = MetricPlumed2([''])  # to test exceptions

    # Simlist
    # datadirs=glob(os.path.join(home(), 'data', 'adaptive', 'data', '*' )
    # fsims=simlist(glob(os.path.join(home(), 'data', 'adaptive', 'data', '*', '/')),
    #              os.path.join(home(), 'data', 'adaptive', 'generators', '1','structure.pdb'))

    import doctest
    doctest.testmod()


    pass
