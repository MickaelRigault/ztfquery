***
skyvision.py documentation
***

`ztfquery.skyvision.py` is a module to access the observing logs of ZTF.

# Requirements 
- You need to know the password to skyvision (internal to the ztf collaboration)

*** 

The `skyvision.CompletedLog` object is able to download and store the logs (`$ZTFDATAPATH/skyvision`) and has plenty of convinient functionalities.

# Single Night Logs

Get the **logs of a given night**, say the 1st of February 2020. Remark that the code will check if this logs is already stored and automatically download it if not.
```python
logs = skyvision.CompletedLog.from_date("2020-02-01")
```

Data are stored as `logs.data`

|     | datetime                | date       |   exptime |   totalexptime |   fid |   field |   pid | ra           | dec       |   totaltime |       obsjd |
|----:|:------------------------|:-----------|----------:|---------------:|------:|--------:|------:|:-------------|:----------|------------:|------------:|
|   1 | 2020-02-01T02:16:18.868 | 2020-02-01 |        30 |        154.423 |     1 |     447 |     1 | +00:20:57.39 | +04:33:00 |     154.423 | 2.45888e+06 |
|   2 | 2020-02-01T02:17:03.249 | 2020-02-01 |        30 |         44.472 |     1 |     603 |     1 | +01:34:9.22  | +26:09:00 |      44.472 | 2.45888e+06 |
|   3 | 2020-02-01T02:17:43.190 | 2020-02-01 |        30 |         39.943 |     1 |     652 |     1 | +02:06:50.27 | +33:21:00 |      39.943 | 2.45888e+06 |
...
| 838 | 2020-02-01T13:47:44.876 | 2020-02-01 |        90 |         98.796 |     3 |     759 |     2 | +14:57:5.74  | +47:45:00 |      98.796 | 2.45888e+06 |
| 839 | 2020-02-01T13:49:26.419 | 2020-02-01 |        90 |        101.633 |     3 |     823 |     2 | +15:33:20    | +62:09:00 |     101.633 | 2.45888e+06 |

To visualize the field observed run:
```python
logs.show_gri_fields(title="2020-02-01")`:
```

<p align="left">
  <img src="images/completedlogs_20200201.png" width="700" title="hover text">
</p>

### Methods

The `CompletedLog` object have several convenient pre-built method that enables you to interact with the dataframe stored as `logs.data`.
For instance:
  - get the entries when a (or list of) field(s)  was observed: `logs.get_when_field_observed(fieldid)`
  - get the entries when a target  was observed: `logs.get_when_target_observed([ra,dec])`
  - get a filtered version of the data: `logs.get_filter(SEE OPTIONS)`, use as `logs.data[logs.get_filter(OPTIONS)]`


# Multiple Day Logs

This works the very same way as for a single day, simply do:
```python
logs = skyvision.CompletedLog.from_date(["2020-02-01","2020-07-03"])
```

# Logs between time range

This works again the very same way as for a single day, but you can simply loads it as:

```python
logs = skyvision.CompletedLog.from_daterange("2020-02-01",end=None)
```

Then for instance:

```python
logs.get_when_field_observed(456)
```

|     | datetime                | date       |   exptime |   fid |   field |   pid | ra           | dec       |   totaltime |       obsjd |
|----:|:------------------------|:-----------|----------:|------:|--------:|------:|:-------------|:----------|------------:|------------:|
| 176 | 2020-02-01T04:18:41.870 | 2020-02-01 |        30 |     2 |     456 |     1 | +04:33:20.89 | +04:33:00 |      38.87  | 2.45888e+06 |
| 106 | 2020-02-05T03:32:29.182 | 2020-02-05 |        30 |     2 |     456 |     1 | +04:33:20.89 | +04:33:00 |      38.793 | 2.45888e+06 |
| 162 | 2020-02-07T04:12:57.286 | 2020-02-07 |        30 |     1 |     456 |     1 | +04:33:20.89 | +04:33:00 |      38.855 | 2.45889e+06 |
...
|  59 | 2020-03-05T03:22:17.415 | 2020-03-05 |        30 |     1 |     456 |     1 | +04:33:20.89 | +04:33:00 |      41.628 | 2.45891e+06 |
| 137 | 2020-03-05T04:15:26.386 | 2020-03-05 |        30 |     2 |     456 |     1 | +04:33:20.89 | +04:33:00 |      39.146 | 2.45891e+06 |
| 122 | 2020-03-09T04:25:58.799 | 2020-03-09 |        30 |     1 |     456 |     1 | +04:33:20.89 | +04:33:00 |      38.855 | 2.45892e+06 |


You can also count the number of time each program when used to observe a given filter (say 'ztf:r' so fid=2):
```python
logs.get_count("pid", fid=2)
pid
1    24986
2    13661
3    14249
```

As a final example of things you can do, watch what ZTF did:

![](images/logsurvey.gif)

# Bulk downloading

First download the logs on your laptop by doing (in this example since the 1st of May 2018):
```python
from ztfquery import skyvision
skyvision.download_timerange_log("2018-05-01", which="completed", nprocess=4)
```
A progress bar is prompted. It roughly takes ten to twenty seconds per year. You will only need to do that once.

Alternatively, the necessary logs are automatically downloaded when loading a `CompletedLog` object if necessary.
