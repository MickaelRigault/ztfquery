# Reference Images

If you want to know if a given field (say 400) already have there reference images use:
```python
from ztfquery import fields
fields.has_field_reference(400)
"""
{'zg': True, 'zi': False, 'zr': True}
"""
```

If you want the list of all field that have, say a I-band image:
```python
from ztfquery import fields
fields.get_fields_with_band_reference("zi")
"""
441,  442,  516,  517,  518,  519,  520,  522,  523,  524,  525,
526,  527,  528,  530,  531,  532,  534,  544,  547,  549,  550,
564,  565,  566,  567,  568,  569,  570,  571,  572,  573,  574,
575,  576,  577,  578,  579,  580,  581,  582,  583,  584,  585,
586,  596,  597,  613,  615,  616,  617,  618,  619,  620,  621,
622,  623,  624,  625,  626,  627,  628,  629,  630,  631,  632,
633,  634,  635,  645,  646,  660,...
"""
```
