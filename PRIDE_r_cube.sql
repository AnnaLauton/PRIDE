create or replace function pride_r(
    p_table             regclass,	-- data table
    p_id_col            text,		-- identifier column
    p_feat_col          text,      -- feature vector (cube) column
    p_pivots_table      regclass,  -- pivot table
    p_pivot_pos_col     text,      -- pivot ordering column
    p_pivot_cube_col    text,      -- pivot feature vector (cube) column
    p_center_id         bigint,    -- query center id
    p_r                 float8,    -- query radius
    p_m                 int        -- rings are indexed from 0 to m
)
returns table (
    id              bigint,			-- returned object identifier
    dist_center     float8,			-- distance from the returned object to the query center
    pivot_radical	int[]			-- pivot-ring radical of the returned object
)
language plpgsql
as $$
declare
    sql text;
begin
    sql := format($fmt$
with

-- Query parameters
params as (
    select
        %1$L::bigint as center_id,
        %2$L::float8 as r,
        %3$L::int    as m
),

-- Query center feature vector
q as (
    select c.%5$I as center
    from %4$s c
    join params on c.%6$I = params.center_id
),

-- Precomputed maximum distance between dataset objects and pivots
metadata as (
    select max_dist
    from dataset_pivot_metadata
    where data_table   = '%4$s'::regclass
      and pivots_table = '%7$s'::regclass
),

-- Ring thickness used to discretize distances to pivots
ring_params as (
    select
        p.center_id,
        p.r,
        p.m,
        md.max_dist,
        (md.max_dist / (p.m + 1))::float8 as band_size
    from params p
    cross join metadata md
),

-- Ordered pivot set
pivots as (
    select
        p.%8$I::int as pivot_pos,
        p.%9$I      as pivot_cube
    from %7$s p
    order by p.%8$I
),

-- Objects within the query radius
candidates as (
    select
        d.%6$I::bigint as id,
        d.%5$I         as feat,
        (q.center <-> d.%5$I)::float8 as dist_center
    from %4$s d
    cross join q
    cross join params
    where d.%6$I <> params.center_id
      and d.%5$I && cube_enlarge(q.center, params.r, cube_dim(q.center))
      and (q.center <-> d.%5$I) <= params.r
),

-- Distances from candidates to pivots
candidate_pivots as (
    select
        c.id,
        c.dist_center,
        p.pivot_pos,
        (c.feat <-> p.pivot_cube)::float8 as dist_cp
    from candidates c
    cross join pivots p
),

-- Ring index of each candidate with respect to each pivot
candidate_rings as (
    select
        cp.id,
        cp.dist_center,
        cp.pivot_pos,
        least(
            floor(cp.dist_cp / nullif(rp.band_size, 0))::int,
            rp.m
        ) as ring_id
    from candidate_pivots cp
    cross join ring_params rp
),

-- Radical: ordered vector of pivot-ring indices
radicals as (
    select
        id,
        dist_center,
        array_agg(ring_id order by pivot_pos) as pivot_radical
    from candidate_rings
    group by id, dist_center
),

-- Keep the closest object to the query center for each radical
ranked as (
    select
        s.id,
        s.dist_center,
        s.pivot_radical,
        row_number() over (
            partition by s.pivot_radical
            order by s.dist_center, s.id
        ) as rn
    from radicals s
)

select
    id,
    dist_center,
    pivot_radical
from ranked
where rn = 1
order by dist_center, id
$fmt$,
        p_center_id,       -- %1$L
        p_r,               -- %2$L
        p_m,               -- %3$L
        p_table,           -- %4$s
        p_feat_col,        -- %5$I
        p_id_col,          -- %6$I
        p_pivots_table,    -- %7$s
        p_pivot_pos_col,   -- %8$I
        p_pivot_cube_col   -- %9$I
    );

    return query execute sql;
end;
$$;
