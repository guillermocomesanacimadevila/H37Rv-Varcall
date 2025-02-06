-- vcfs_1
-- vcfs_2
-- vcfs_3

SELECT COUNT(*) FROM vcfs_1; -- 1840 entries

SELECT COUNT(*) FROM vcfs_2; -- 1995 entries

SELECT COUNT(*) FROM vcfs_3; -- 2493 entries

SELECT * FROM vcfs_1;

-- Let´s have a look at if there are any INDELs in samples 1-3

-- Sample 1 - YES! Several INDELs on the Reference allele

SELECT REF, COUNT(*) AS REF_Count
FROM vcfs_1
GROUP BY REF
ORDER BY REF_COUNT;

SELECT ALT, COUNT(*) AS ALT_Count
FROM vcfs_1
GROUP BY ALT
ORDER BY ALT_Count;

-- Sample 2 -- Just all 4 DNA nucleotides

SELECT REF, COUNT(*) AS REF_Count
FROM vcfs_2
GROUP BY REF
ORDER BY REF_COUNT;

SELECT ALT, COUNT(*) AS ALT_Count
FROM vcfs_2
GROUP BY ALT
ORDER BY ALT_Count;

-- Sample 3 -- Once again, just ACTG

SELECT REF, COUNT(*) AS REF_Count
FROM vcfs_3
GROUP BY REF
ORDER BY REF_COUNT;

SELECT ALT, COUNT(*) AS ALT_Count
FROM vcfs_3
GROUP BY ALT
ORDER BY ALT_Count;

/** So we know sample 1 needs cleaning and it contains INDELs **/

SELECT COUNT(*)
FROM vcfs_1
WHERE LENGTH(ALT) > 1;

-- 42, but this is not accounting for all insertions and deletions

SELECT 
    COUNT(*) AS total_variants,
    COUNT(CASE WHEN LENGTH(ALT) > LENGTH(REF) THEN 1 END) AS insertions,
    COUNT(CASE WHEN LENGTH(ALT) < LENGTH(REF) THEN 1 END) AS deletions
FROM vcfs_1; -- 57 INDELs

-- Clean VCF

SELECT *
FROM vcfs_1
WHERE LENGTH(ALT) <> LENGTH(REF);