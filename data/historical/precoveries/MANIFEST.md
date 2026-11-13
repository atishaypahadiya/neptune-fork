# Held-out precoveries of Neptune (M6e)

Held out: no fitting module may read this directory. `tests/test_a6d.py`
(`test_no_reference_to_a_precovery_file`) fails if a module on the fit path names it.
Sourced 12 November 2026. Nothing has been compared with any posterior.

## Lalande, 8 and 10 May 1795: sourced

V. Mauvais, "Sur une observation inédite de la nouvelle planète", *Comptes rendus* 24 (1847),
666–669. The reduced apparent places are in the table on p. 668; the observation times were
corrected by an erratum in the same volume. A second printing: Mauvais's letter of 5 July
1847, *Astronomische Nachrichten* 26, No. 607, col. 97.

| quantity | 8 May 1795 | 10 May 1795 |
|---|---|---|
| Paris mean time (erratum) | 11h 5m 11s | 10h 57m 9s |
| apparent right ascension | 212° 59′ 35″.0 | 212° 56′ 36″.3 |
| apparent declination | −11° 20′ 39″.1 | −11° 19′ 38″.8 |

Apparent places: equinox and equator of date, aberration and nutation not removed. Mauvais
states no uncertainty. He says the places "ne contiennent plus d'autres erreurs que celles
qui sont inhérentes à l'observation elle-même". The raw transits and zenith distances are on
p. 667 and were not transcribed. The printed differences (−2′ 58″.7, +1′ 0″.3) equal the
differences of the rows.

Two entries were made from three page bands, the second from the bands and a column
specification alone: 9 rows, 117 cells, 0 disagreements. The two printings agree in both
places and in the time of 8 May. They differ in the time of 10 May (57m against 59m); the
erratum's value is adopted and the reason is in the misprint log.

| file | content | SHA-256 |
|---|---|---|
| `lalande_1795.csv` | the two adopted places | `22484a07c4e8ee568ca3b84d4f35b16a9006aeb16d935a65f5ce0b095885849a` |
| `lalande_1795_misprints.csv` | four entries: two errata applied, two logged | `6d5492cc56a51ec4a1d607c05ffbef0ff00c0b828cb2def0f7508ddcc852a351` |
| `passes/lalande_1795_pass1.csv` | first entry | `d07bb8ce0969411a2251f162204308a594b02629414a702b5c625be3f7504304` |
| `passes/lalande_1795_pass2.csv` | second entry | `f9e3fb1fd18796bf06b27b9fa768ebdc347c15501995d2b51914b24e6dad3da6` |

Sources (cache, not committed; fetched 12 November 2026):

| file | source | SHA-256 |
|---|---|---|
| `data/cache/adams_challis/prediscovery/cr24_leaf_n676_w2000.jpg` | archive.org `comptesrendusheb24acad`, leaf n676 (p. 667) | `08e55ad2b20c77798e976a1f8f2cb773a28b568073abfa1d6803bb151daed3b8` |
| `.../cr24_leaf_n677_w2000.jpg` | leaf n677 (p. 668) | `2ae8a0f70eac67193ebc856b1b13d4e393f69b880e3910e38cca7ef4d57abd3b` |
| `.../cr24_leaf_n877_w2000.jpg` | leaf n877 (the errata; its page number was not read) | `87d4cdf75fd4ea4fae6650b18be1c289e1a5d12172177befc787f3ef4dd4bba0` |
| `.../1847AN.....26...97M.pdf` | ADS scan of *Astronomische Nachrichten* 26, col. 97 | `b00a0eae7f9e21aface6d5872a05c458ba0eaaf7f98caada682f35902465ea2f` |
| `data/cache/precoveries/bands/cr24_p668_table.png` | band cut from leaf n677 | `fab2b3fb0ff5eb0251e53ce89e8bd78a8c89cff5a52760aafa383ff15e7b09a5` |
| `.../cr24_errata.png` | band cut from leaf n877 | `70ec7adcc5f4f3c7e7fc0037a80fba049db653ad6b55efbbc1b72b1bb234d21b` |
| `.../an26_col97.png` | band cut from the ADS scan | `1329c27fc4d25ac4009f7df5cb3d3879cbc42e491c13ba3753a5e5b0a3112f36` |

## Galileo, 28 December 1612 and 27 and 28 January 1613: not sourced

Kowal and Drake 1980, *Nature* 287, 311 has not been read beyond its abstract (the journal
needs a subscription or a library). The ADS scan of Standish and Nobili 1997, *Baltic
Astronomy* 6, 97 serves its first page only by script. What is known from Standish's summary
(*Highlights of Astronomy* 12, 326–327): the record is a drawing, not a measured place. On
28 January 1613 Neptune is drawn on the line from Jupiter through a catalogued star, and its
drawn separation from that star disagrees with modern ephemerides by an amount "corresponding
to about 53 seconds of arc in Neptune's orbital longitude"; Standish concludes that "at
best" it "can be used as a one-dimensional observation". No file is committed for Galileo.
Needed: the two papers, read in a browser or a library, for the drawn offsets in Jovian radii.
