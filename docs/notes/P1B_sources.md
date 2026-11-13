# Phase 1b, step 0: source inventory for the Adams–Challis study

Date: 2026-11-10. Scope: locate and verify the sources; no analysis, no fitting. Scans are cached
in `data/cache/adams_challis/` (git-ignored); their URLs, licences and SHA-256 hashes are in
`data/historical/adams_challis/MANIFEST.md`, which also says why only one scan (Challis's MNRAS
account, CC0 on Zenodo) is committed.

Every fact below carries one of three marks:

- **[P]** primary-confirmed: read in a scan of the 1846–47 printing (or of the 1904 facsimile
  description) named in the entry;
- **[S]** secondary-only: stated by a later author and not yet seen in a primary source;
- **[U]** unverified: claimed somewhere, not found in anything examined here.

Primary texts were read from text layers or OCR (pypdf, archive.org OCR); numbers quoted from OCR
were checked against a second printing where one exists (MNRAS against Memoirs, Memoirs against the
Nautical Almanac pamphlet). Where the printings differ, both readings are given. No number here has
yet been checked against a page image except the page numbers used for the page mapping.

## 1. Challis's account of the search

**Citation.** J. Challis, "Account of Observations at the Cambridge Observatory for detecting the
Planet exterior to Uranus", *Monthly Notices of the Royal Astronomical Society* 7 (9), 145–149,
13 November 1846, DOI 10.1093/mnras/7.9.145. Full version: "An Account of Observations undertaken in
search of the Planet discovered at Berlin on Sept. 23, 1846", *Memoirs of the Royal Astronomical
Society* 16 (1847), 415–426, read 13 November 1846, dated Cambridge Observatory 12 November 1846.

**Scans.** MNRAS: Zenodo 1431903 (CC0, committed copy), ADS 1846MNRAS...7..145C (text layer).
Memoirs: ADS 1847MmRAS..16..415C (12 image pages, no text), and the archive.org Google scan of the
whole volume (`memoirsroyalast01socigoog`, printed 415 = PDF page 450, checked on the image; OCR in
the `_djvu.txt`). Licences: MANIFEST.

**Relation between the two printings [P].** The MNRAS text is the Society's abstract in the third
person ("Professor Challis says"); the Memoirs text is Challis's own first-person narrative with
more detail (the plan paper, the two methods, assistants, the Sept. 29 position, the comparison-star
tables). The study should cite the Memoirs and use MNRAS as the second witness.

**What it contains that the study needs.**

Search plan (Airy's, as Challis reports it) [P, both printings]:
- Airy's suggestion to search came in a note dated 9 July 1846; the plan paper "Suggestions for the
  Examination of a portion of the Heavens…" followed shortly after (Airy's own account dates its
  despatch 13 July; see § 3).
- Zone: centred on the ecliptic at longitude 325°, "the point marked out by Mr. Adams and
  M. Le Verrier", 15° of longitude each way and 5° of latitude each way; divided into declination
  intervals 1° broad with the RA extent marked for each [P, Memoirs 416–417].
- Power not lower than 120, with which a zone of 15′ would be well seen; transits at one wire;
  NPD estimated against horizontal wires 3′ apart; whole tract swept three times, each sweep
  finished before the next; magnitudes always noted; results catalogued and charted.
- Estimate: 80 zones averaging one hour, three sweeps, about 300 hours of observing.

Challis's modifications [P]:
- Power 166 (the lowest the micrometer eyepiece took), field "not more than 9′"; transits at the
  toothed edge of the micrometer comb; declinations from the comb to a fourth (sometimes a tenth) of
  a revolution, one revolution = 17″ (Memoirs 417–418; MNRAS 145 agrees). Challis's letter of
  12 October 1846 printed in Airy's account says "a power of 160" [P, MNRAS 7, 143]: the two
  printed figures differ.
- Limit: every star visible with the field illuminated enough to show the comb, i.e. to the 11th
  magnitude inclusive [P]. Adams had told him in conversation that a mass about three times that of
  Uranus implied brightness not below the 9th magnitude; Challis kept to the 11th [P, Memoirs 423].
- Two methods, "telescope fixed" (one assistant with memorandum book and chronometer) and
  "telescope moving" (hour circle driven by clock, two assistants, Mr Morgan and Mr Breen)
  [P, Memoirs 418–419]. Plan: go over the tract first fixed, then moving; a third pass only if needed.
- On moonlit nights only "telescope fixed" was used. When a cloud interrupted, a line was drawn in
  the memorandum book and the clock set going; after the cloud the last star was re-observed
  [P, Memoirs 420]. Reference stars taken in each series with the 80′ sector [P].

Berlin charts [P, Memoirs 421]: "In a note, dated July 21, Mr. Airy remarked that but a small portion
of the Berlin maps applied to this inquiry. I found accordingly, in the library of the University of
Cambridge, hour xxii., which extends over a small part of the space intended to be scrutinised; but I
did not find hour xxi., and consequently concluded that it was not published. If I had had this map,
a first sweep would have been unnecessary." MNRAS 146 has the shorter form ("I was ignorant of its
existence, as that hour was not to be found in the University Library"). So the fact that Challis had
access to Hour XXII in the University Library is primary; that he "owned" it [S, Kollerstrom] is not.

Sweep log, as far as it is in print [P, Memoirs 421–423; MNRAS 146–147; Challis to Airy 12 Oct 1846
in MNRAS 7, 142–143]:

| date (1846) | method | what was done |
|---|---|---|
| 29 July | mostly fixed, trials | telescope on the most probable place from Adams's ephemeris (325°) |
| 30 July | moving | zone 9′ broad, all stars to the 11th magnitude |
| 4 August | sector moved in declination | brighter stars over about 70′ (Memoirs: 80′ range of the sector) as reference points; same declination-rod reading as 29–30 July; the planet among them |
| 12 August | fixed | the zone of 30 July repeated, longer in RA; the planet is No. 49 of this series |
| "probably the next day" | — | partial comparison of 12 Aug against 30 July, stopped at No. 39 at a cloud line in the memorandum book; No. 49 (8th magnitude) was the planet |
| 12 Aug – 28 Sept | principally fixed, part repeated moving | "observed subsequently earlier in right ascension"; "continued during two months at all available opportunities"; "some thousands of stars"; did not meet the planet again, "as I took positions too early in right ascension" |
| 29 September | — | after first reading Le Verrier's paper of 31 August that day: a broad sweep in declination within Le Verrier's limits, 3½ hours, nearly 300 stars, one noted "It seems to have a disk" (the planet) |
| 30 September | — | no observations (moon) |
| 1 October | — | news of Galle's discovery arrives |

- Totals [P]: "above 3000" positions recorded (Memoirs 423; MNRAS 147 "3000 positions of stars"),
  contained "in but a small portion of the space intended to be explored". Preparations for mapping
  the stars were made but not resumed.
- Not in print: a night-by-night list of dates, the RA and declination limits of each zone, and the
  star count per night. Only the dates above are given. The actual area covered, needed for the
  "would Challis have found it" question, is in the observing (memorandum) books only (§ 4).
- The MNRAS abstract says the comparison was of "July 30 and August 12" [P]; Krajnović 2016 says
  "between the nights of 4 and 12 August" [S, contradicted by the primary].

Positions [P, Memoirs 424–426; MNRAS 148]: Aug. 4 13h 36m 25s GMT, RA 21h 58m 14.70s, NPD
102° 57′ 32.2″; Aug. 12 13h 3m 26s, RA 21h 57m 26.13s, NPD 103° 2′ 0.2″ (means over the comparison
stars tabulated in the Memoirs; NPD errors about 3″ and 4″); Sept. 29 10h 10m, RA 21h 52m 48.9s,
NPD 103° 26′ 30″ from one star, B.A.C. 7599 (Memoirs only). Adams's elements from these and
observations to 13 October: distance 30.05, heliocentric longitude 326° 39′ on 4 August, descending
node 309° 43′, inclination 1° 45′ [P].

Adams's ephemeris for Challis [P, Memoirs 421–422; MNRAS 146–147]: RA and declination for every
twentieth day from 20 July to 8 October and for every fifth degree of heliocentric longitude from
315° to 335°, planet assumed in the ecliptic; the 325° column is printed: 20 July 21h 51.3m −13° 0′;
9 Aug 21h 49.5m −13° 9′; 29 Aug 21h 47.9m −13° 18′; 18 Sept 21h 46.5m −13° 26′; 8 Oct 21h 45.1m
−13° 33′. The other four columns are not printed; the manuscript is the "July ephemeris" (§ 4, § 5).

## 2. Adams's memoir

**Citations.**
- Abstract: J. C. Adams, "An Explanation of the observed Irregularities in the Motion of Uranus, on
  the Hypothesis of Disturbance caused by a more distant Planet; with a Determination of the Mass,
  Orbit, and Position of the disturbing Body", MNRAS 7 (9), 149–152, 13 November 1846,
  DOI 10.1093/mnras/7.9.149 [P, ADS 1846MNRAS...7..149A, text layer].
- Full memoir, first printing: Appendix to *The Nautical Almanac and Astronomical Ephemeris for the
  year 1851* (London, 1847), "On the Perturbations of Uranus", pp. 265–293 [S for this pagination:
  Kollerstrom 2006 and library records; the HathiTrust catalogue record 000552896 lists 1851 volumes
  but refused automated access, so the volume scan was not opened]. The separate printing (W. Clowes
  & Sons, 1846, "From the Appendix to the Nautical Almanac for the year 1851", 31 pp. with its own
  pagination, §§ 1–61, signed St John's College 12 November 1846) is on archive.org,
  `explanationobse00Adam` [P]. The French translation's *avertissement* explains the venue: the
  Society's printer was occupied with a long paper on the longitude of Valentia, so the memoir went
  to the Almanac appendix first [P, JMPA 1876, 5].
- Memoirs RAS 16 (1847), 427–460, read 13 November 1846 [P: start page 427 on the image; 34 ADS
  pages give 427–460]. ADS 1847MmRAS..16..427A (no text); archive.org volume, printed 427 = PDF 464.
- French translation: "Explication des irrégularités observées dans le mouvement d'Uranus…",
  *Journal de Mathématiques Pures et Appliquées* 3e sér. 2 (1876), 5–32 and 69–86 (*suite*)
  [P, Numdam `JMPA_1876_3_2__5_0`, `JMPA_1876_3_2__69_0`, text layer]. The translator is not named
  on the first pages read; the heading calls it "Extrait de l'Appendice du Nautical Almanac pour
  l'année 1851".

**Scans.** As above; licences in MANIFEST.

**What it contains that the study needs.**

Chronology of Adams's work, in his own words [P, pamphlet §§ 3–5; Memoirs 428–429; MNRAS 149–150]:
- Memorandum dated 3 July 1841: resolve to investigate after taking his degree.
- 1843: first solution, circular orbit at twice Uranus's mean distance, modern observations only.
- February 1844: Airy sends the reduced Greenwich observations of Uranus.
- Second solution with first-order eccentricity, same distance; "after obtaining several solutions
  differing little from each other" the final values were given to Challis in September 1845 and,
  "slightly corrected", to Airy the following month (Hypothesis I).
- Repeated with mean distance "about 1/30th part less" (Hypothesis II), communicated to Airy "in the
  beginning of September" 1846 (the letter of 2 September), with the inference that the distance
  should be reduced further.
- Orbit from Challis's August positions plus observations to 13 October (in Challis's paper, § 1).

Elements [P]:

| | Hyp. I (a′/a = 2) | Hyp. II (a′/a = 1.942) |
|---|---|---|
| mean longitude at epoch 1810.328 | 269° 25′ | 264° 50′ |
| mean longitude 6 October 1846 | 325° 7′ | 323° 2′ |
| longitude of perihelion (1846) | 315° 57′ | 299° 11′ |
| eccentricity | 0.16103 | 0.120615 |
| mass (Sun = 1) | 0.0001656 | 0.00015003 |

(pamphlet §§ 31, 48; MNRAS 151 agrees; the 2 September letter in Airy's account gives the Hyp. I
mean longitude for 1 October 1846 as 325° 8′ and Hyp. II as 323° 2′, with masses 0.00016563 and
0.00015003.) The October 1845 paper left at Greenwich gave: mean distance 38.4 ("assumed nearly in
accordance with Bode's law"), mean sidereal motion 1° 30′.9 per 365.25 days, mean longitude
1 October 1845 323° 34′, perihelion 315° 55′, eccentricity 0.1610, mass 0.0001656 [P, Airy's
account, letter No. 11, MNRAS 7, 129].

Residuals [P, pamphlet §§ 52–53; MNRAS 151–152]: maximum residuals 1712 +6.7″/+6.3″, 1715 −6.8/−6.6,
1753 +5.7/+5.2, 1771 +11.8/+12.8 (single observation); after 1771 to 1840 the largest is 2.35″;
1843 +7.11/+5.77, 1844 +8.79/+7.05, 1845 +12.40/+10.18 (Hyp. I / Hyp. II, observation minus theory,
mean longitude). 1690 (Flamsteed, excluded from the equations): +44.5″ / +50.0″. Radius-vector
corrections 1834/1840/1846: +0.00505/+0.00722/+0.00868 (I), +0.00492/+0.00696/+0.00825 (II).

The 2 September 1846 inference [P, Airy's account, letter No. 20, MNRAS 7, 137–139; pamphlet § 53]:
from the 1843 errors (6.84″ and 5.50″ in the letter; 5 : 4 over 1843–45 in the memoir) Adams
inferred a/a′ = 0.57 (sin 35° = 0.574) "and the corresponding mean longitude on the 1st October,
1846, would be about 315° 20′, which I am inclined to think is not far from the truth". The memoir
(§ 54) withdraws the direction of this extrapolation: the longitude at a smaller distance "would
probably soon arrive at a minimum, and afterwards begin again to increase", and he "hastily inferred"
otherwise in the letter. Kollerstrom (note 129) cites Adams to Airy, 15 October 1846, to the same
effect [S]. Neptune's true heliocentric longitude on 23 September 1846 was 326° 57′ [S, Kollerstrom
2006 App. II], so 315° 20′ is about 11½° off; that figure is derived, not quoted.

"Six solutions" [U]: no primary source counts the solutions. The primary record gives: 1843; the
"several solutions differing little from each other" of 1844–45 (no count); Sept. 1845 (Challis) and
Oct. 1845 (Airy, "slightly corrected"), Hyp. I; Hyp. II of 2 Sept. 1846; the 0.574 extrapolation of
the same letter; the November 1846 orbit from Challis's observations. The dioi.org "Adams Dated
Computations" page lists two dated solutions in the manuscripts, 18 September 1845 (−50° 34′, the
difference of mean longitudes at the 1810 epoch) and 20 August 1846 (Hyp. II, −42° 52′), plus dated
radius-vector work (12 Nov., 28 Nov., 15–16 Dec., 24 Dec. 1845; 1 Sept. 1846) [S]. The published
inventory of the manuscripts is Sampson 1904 (§ 5); reading it is the step-1 route to a complete
list.

## 3. Airy's account

**Citation.** G. B. Airy, "Account of some Circumstances historically connected with the Discovery of
the Planet exterior to Uranus", read 13 November 1846; MNRAS 7 (9), 121–144 [P, ADS
1846MNRAS...7..121A, text layer; this ADS item is the whole report of the 13 November meeting, and
Crossref lists MNRAS 7, 121 under another title, so cite by page]; Memoirs RAS 16 (1847), 385–414
[P: 385 on the image; 30 ADS pages]; reprinted in *Astronomische Nachrichten* 25, no. 585 (1847),
DOI 10.1002/asna.18470251002 [P, Zenodo 1424621, CC0, 14 pp., text layer].

**What it contains that the study needs** [P, MNRAS pagination]:
- Numbered letters. No. 10, Airy to Challis 29 Sept. 1845 (p. 129); No. 11, Adams's October 1845
  paper with the elements above (129–130) and its table of residuals; Airy to Adams 5 Nov. 1845 (the
  radius-vector question, 130); Airy to Le Verrier 26 June 1846 (132); No. 15, Airy to Challis from
  the Deanery, Ely, 9 July 1846, proposing the Northumberland telescope and an assistant (135);
  No. 16, 13 July 1846, sending the "Suggestions", which "contemplated the examination of a part of
  the heavens 30° long, in the direction of the ecliptic, and 10° broad" (136); No. 17, Challis to
  Airy 18 July ("I purpose to carry the sweep to the extent you recommend"); No. 18, Challis to Main
  7 Aug.; No. 19, Challis to Airy 2 Sept. (slow progress, stars to magnitude 10–11, "many more
  observations than I can take this year", 137); No. 20, Adams to Airy 2 Sept. 1846 (137–139);
  No. 23, Challis to Airy 12 Oct. 1846, the narrative quoted in § 1 (142–143).
- Le Verrier's 31 August elements as Airy summarises them (140–141): semi-axis 36.154, eccentricity
  0.10761, perihelion 284° 45′, mean longitude 1 Jan. 1847 318° 47′, mass 1/9300, true heliocentric
  longitude 1 Jan. 1847 326° 32′, distance 33.06.
- The text of the "Suggestions" paper is not printed; Challis paraphrases it (§ 1). The original is
  in the archives (§ 4). Airy's notes of 12 and 21 July 1846 to Challis about the Berlin maps are
  known only from Challis's Memoirs sentence and from Kollerstrom's citation "Airy to Challis,
  12 July 1846, COA 4 (attachment p. 4)" [S].
- "Airy steered the search toward Le Verrier's position" [S]: what is primary is that the zone centre
  was 325°, described by Challis as marked out by both Adams and Le Verrier, that Airy's letter of
  9 July followed Le Verrier's letter of 28 June by "about a week", and that Challis changed his
  method on 29 September on reading Le Verrier's 31 August paper. The stronger reading is
  Kollerstrom's and Rawlins's.

## 4. Archives

Everything in this section is catalogue-confirmed only (the catalogue entry was read; nothing
inside the files has been seen), unless marked otherwise.

**Cambridge University Library, University Archives (Cambridge Observatory).**
- *Discovery of Neptune papers, 1845–1976*, GBR/0265/UA/Obsy D.3 (ArchiveSearch
  `archival_objects/23027`): "correspondence and observations; photocopy of article on G.B. Airy in
  the 'Eagle', Easter 1976"; extent "4 folders of letters and 2 notebooks"; "List of letters compiled
  by D.W. Dewhirst housed with originals"; open to holders of a Manuscripts Reading Room ticket.
  Whether the two notebooks are the 1846 sweep memorandum books is not stated [U]. Kollerstrom cites
  Cambridge Observatory material as "COA n" [S].
- The Observatory fonds (UA/Obsy, `archival_objects/8787`) has an "Observations and calculations
  records" series (`archival_objects/10104`); the records pages returned empty bodies to automated
  requests, so the sub-series listing (where the Northumberland equatorial observing books of 1846
  would sit) is unread [U]. A reading-room enquiry is needed.
- *Observatory, 1818–1878*, UA/CUR 29 (`archival_objects/16087`): Syndicate papers; Challis's report
  to the Observatory Syndicate of 12 December 1846 [S, dioi.org "Challis' Unseen Testimony"] would be
  here or in UA/Obsy.

**Cambridge University Library, Royal Greenwich Observatory Archives (Airy papers).**
- *Papers on Neptune and the minor planets, 1834–1884*, GBR/0180/RGO 6/96a
  (`archival_objects/98287`): "Correspondence of G.B. Airy (including many copy letters) and related
  papers concerning planets, including circulars, ephemerides, observations, announcements of
  discoveries, discussion of names…"; "The volume is currently undergoing conservation work and is
  unavailable to readers, although a microfilm copy may be consulted". The microfilm is RGO 36/2/70,
  "Microfilm of RGO 6/96–98", made by the Public Record Office in 1979. This is the "Neptune file":
  Airy's 1846 letters, the "Suggestions" paper, Challis's letters of 18 July, 7 Aug., 2 Sept. and
  12 Oct., Adams's letters of 2 Sept. and 15 Oct., and (per Kollerstrom, note 119, "mis-filed")
  Adams's July 1846 ephemeris [S for the contents].
- RGO 6/97 (1848–50), RGO 7/247 (*Miscellaneous correspondence, 1846–1897*, `archival_objects/99289`:
  "duplicate copies of collected correspondence on the investigations of J.C. Adams and
  U.J.J. Le Verrier into the existence of Neptune and in particular of the ensuing controversy over
  the name of the planet") and RGO 7/248 are the duplicates.
- Cambridge Digital Library has digitised parts of the Airy papers (RGO 6); whether RGO 6/96a is
  among them was not established [U].

**St John's College, Cambridge.**
- "The 'Neptune Notebook' of John Couch Adams, 1840s (MS W.16)" is listed among the Old Library's
  highlights; no digitisation statement on the page [P for the listing].
- Papers of John Couch Adams (Archives Hub 736b6931-3eea-3de3-b245-eafbd807ad5d; boxes 1, 20,
  43–44 relevant per the earlier reading of the hub record, which refused automated access today)
  [S]. Sampson 1904 (§ 5) describes "a volume of manuscripts … presented by Mrs. Adams to the
  Library of St. John's College" containing "all, or very nearly all" that Adams wrote on the
  perturbations of Uranus from the 1841 memorandum to the discovery [P, Sampson 143], and his
  catalogue uses the sigla E III, E VII, B VIII [S, via Kollerstrom's notes].
- MacAlister's transcripts of the correspondence ("McA box:folder" in Kollerstrom) [S].

**Berlin chart, Hora XXII.** Argelander's sheet of the Akademische Sternkarten, completed 1832,
published 1833 [S, Krajnović 2021 Fig. 6.10 caption; Kollerstrom 2006, 26]. A copy is in the library
of the Leibniz-Institut für Astrophysik Potsdam, which supplied the figure [S]; copies existed in
Cambridge and Paris in 1846 [P for Cambridge, § 1]. Coverage [S]: each sheet is one hour of RA wide
with a 4-minute (1°) overlap on each side (Bessel's 1825 proposal as quoted in the 1859 catalogue),
declination ±15°, epoch 1800; Hour XXII therefore reaches to 21h 56m, and Neptune in August 1846 was
"passing through the rightmost quadrants (between 21h 56m and 22h) at the declination of
approximately −13°". Against the primary positions: Challis's 4 and 12 August places (21h 58m 15s and
21h 57m 26s) fall on the sheet; Adams's predicted 325° place (21h 49.5m on 9 August) does not. No
scan of the sheet has been obtained [U]; the chart's actual limits need the sheet or the 1859
catalogue.

## 5. Secondary sources (cross-checking only)

- D. Krajnović, "The contrivance of Neptune", *Astronomy & Geophysics* 57 (5), 5.28–5.34 (2016),
  DOI 10.1093/astrogeo/atw183; arXiv:1610.06424v2 (cached). States: search started 29 July; Challis
  recorded the planet twice before and once after the discovery; "comparing observations between the
  nights of 4 and 12 August, he stopped at star number 39" (the primary says 30 July and 12 August);
  "it is remarkable that Challis had Hora XXII chart, adjacent and partially overlapping the map in
  Berlin"; Adams's new orbit from Challis's three observations. No count of Adams's solutions and no
  315⅓° figure in this text.
- D. Krajnović, "'That star is not on the map': the German side of the discovery", chapter in
  *Neptune: From Grand Discovery to a World Revealed* (Springer 2021, pp. 185–244, doi:10.1007/978-3-030-54218-4_6), arXiv:2108.06305 (CC BY 4.0, cached). Source for the
  Hora XXII statements above and for the Berlin charts' design (24 one-hour sheets, ±15°, 4-minute
  overlaps, grids of 510 one-degree squares).
- N. Kollerstrom, "An hiatus in history: the British claim for Neptune's co-prediction, 1845–1846",
  *History of Science* 44 (1), 1–28 (Part 1, cached from a mirror) and 44 (3), 349–371 (Part 2,
  DOI 10.1177/007327530604400303, paywalled, not obtained). Part 1 gives: the sweep "began at Cambridge
  on 29 July 1846 … and it lasted two months"; Challis's 17 October 1846 *Athenaeum* letter (first
  public account of the search method; not obtained); Challis saw Neptune on 4 and 12 August; "In
  Cambridge, Challis had the older 'Hour 22' Berlin map published in 1833, produced by Argelander";
  Neptune moved from Hour 22 onto Hour 21 "towards the end of August", so four weeks of the search
  were on the map Challis had; Airy's 12 July letter alluded to the Hour 22 map (note 184, COA 4);
  Appendix I (distances: Le Verrier 36.2, Adams Hyp. I 38.4, Hyp. II 37.5 AU) and Appendix II
  (longitudes converted to 23 September 1846: Adams I 328° 41′, Adams II 329° 27′, Le Verrier June
  324° 35′, August 325° 58′; Neptune 326° 57′); the 2 September 315° 20′ passage with "which I am
  inclined to think is not far from the truth"; Adams's dated manuscript solution of 18 September
  1845 (−50° 34′) from Sampson; note 80: the eccentricity, apse and mean longitude first appear in
  Adams MS E VII, p. 15, "enabling the solution to be sent to Airy on 2 September"; note 33: the
  memorandum-book entry reads "last one seemed to have a disc", part scratched out.
- dioi.org, "Neptune's Discovery: The British Case for Co-Prediction" (N. Kollerstrom), pages under
  `https://dioi.org/kn/neptune/`: `chron.htm` (chronology: 9 July Airy urges the search and Adams
  prepares the ephemeris; 29 July start; 30 Sept. Challis hears the news; 13 Nov. RAS meeting),
  `dated.htm` (the dated computations listed in § 2), `july.htm` (the July ephemeris: five columns
  over 20° at two-week steps 20 July–8 October, lower section circular ecliptic orbit, upper section a
  node and inclination from Wartmann's star, "q = 336°.4 nearly"), `mapl.htm` (Hour 22 to 21h 56m;
  the 17 October positions RA 21h 58m and 21h 57m), `unseen.htm`, `corr.htm`, `sources.htm` (images
  only). The site blocks command-line clients; pages were read through a fetch service. All [S].
- W. Sheehan, N. Kollerstrom, C. B. Waff, "The case of the pilfered planet", *Scientific American*
  291 (6), 92–99 (December 2004), DOI 10.1038/scientificamerican1204-92. Paywalled; not obtained.
- C. J. Brookes, "On the prediction of Neptune", *Celestial Mechanics* 3 (1), 67–80 (1970),
  DOI 10.1007/BF01230433. Paywalled; not obtained. Kollerstrom (note 128) cites p. 78 for how the
  solutions behave near the 5:2 resonance distance.
- R. A. Sampson, "A Description of Adams's Manuscripts on the Perturbations of Uranus", *Memoirs
  RAS* 54 (1904), 143 ff., received 13 December 1901, "with twenty pages of reproductions of MS."
  [P: title page and p. 144 read on the image; ADS 1904MmRAS..54..143S, 48 image pages, cached, no
  text layer]. This is the published catalogue of the St John's volume and the route to every dated
  solution; it is primary for the facsimile plates and the catalogue, secondary for Sampson's
  narrative. Not yet read beyond p. 144.
- D. Rawlins, "British Neptune-disaster file recovered", *DIO* 9.1 (1999) [S, cited by Kollerstrom;
  not obtained; DIO volumes are normally free on dioi.org].

## 6. What is available, what is missing, what needs an archive request

**Available now (cached, hashes in the MANIFEST):** all three accounts of 13 November 1846 in both
MNRAS and Memoirs printings; Airy's AN reprint; Adams's memoir in the separate printing, the Memoirs
and the French translation; Sampson's 1904 description (images); Krajnović 2016 and 2021; Kollerstrom
2006 Part 1; the dioi.org pages (read, not cached); the Cambridge catalogue entries for UA/Obsy D.3,
RGO 6/96a, RGO 36/2/70 and RGO 7/247.

**Missing but obtainable without an archive visit:**
- The *Nautical Almanac for 1851* appendix itself (HathiTrust record 000552896, full-view volumes
  from Michigan and Princeton per the catalogue summary; the pagination 265–293 is unconfirmed).
- Challis's letter in *The Athenaeum*, 17 October 1846 (his first public account; Google Books or
  HathiTrust).
- Kollerstrom 2006 Part 2, Sheehan et al. 2004, Brookes 1970 (library access).
- The 1859 catalogue of the Akademische Sternkarten (sheet limits) and any digitised Hora XXII.
- Sampson 1904 read in full (cached; needs page-by-page reading or OCR).

**Needs an archive request:**
1. Cambridge University Library, Manuscripts Reading Room: UA/Obsy D.3 (the Dewhirst list and the two
   notebooks), and the Observatory's observing books for the Northumberland equatorial, July–October
   1846. Wanted: the night-by-night sweep log (dates, RA and declination limits of each zone, stars
   per night, method), the memorandum-book pages of 30 July and 12 August with the line after No. 39,
   and the 29 September entry "seems to have a disc". This is the one source for the coverage
   question; nothing in print replaces it.
2. Cambridge University Library, RGO 6/96a via microfilm RGO 36/2/70: Airy's "Suggestions" of 13 July
   1846 with the zone table, Airy's notes of 12 and 21 July, Adams's July 1846 ephemeris (all five
   columns), and the letters listed in § 3. Ask whether RGO 6/96a has been digitised since the
   conservation note was written.
3. St John's College Library: MS W.16 and the volume of Adams's Uranus manuscripts described by
   Sampson (E III, E VII, B VIII), for the dated solutions; ask whether a digitised copy exists.
4. Leibniz-Institut für Astrophysik Potsdam library (or Cambridge University Library's copy): a scan
   of Argelander's Hora XXII sheet, for the chart's limits and the stars plotted near 21h 57m, −13°.

## Addendum, step 1 (2026-11-10)

- The "six solutions" marked [U] in § 2 are primary via Sampson (1904), Memoirs RAS 54, 143 ff.,
  whose Table (p. 151) lists six determinations at epoch 1810.328 with dates for five of them. They
  are transcribed in `data/historical/adams_challis/adams_solutions.csv`; see
  `P1B_adams_solutions.md`.
- The 2 September 1846 figure of 315° 20′ is primary but is an extrapolation, not a solution; the
  verdict on the "~315⅓°, ~12° from Neptune" claim is in `P1B_adams_solutions.md`.

## Addendum, step 2 (2026-11-10)

- Correction to § 4, "Berlin chart, Hora XXII": the remark that Challis's 4 and 12 August places
  (21h 58m 15s and 21h 57m 26s) "fall on the sheet" compared right ascensions of date with limits
  stated for the equinox of 1800. In the equinox of 1800 those places are 21h 55m 43s and
  21h 54m 54s, both west of the 21h 56m edge of Hora XXII (with the secondary 4-minute overlap); the
  planet left the sheet on 2–3 August 1846. See `P1B_challis_sweeps.md`, "The Berlin charts". The
  limits remain [S].
- The search log is reconstructed in `data/historical/adams_challis/challis_sweeps.csv`; what the
  printed record cannot give (zone boundaries, the dates between 12 August and 28 September) is listed
  in `P1B_challis_sweeps.md`.
