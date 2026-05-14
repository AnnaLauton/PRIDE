create or replace function pride_k_omni_array(
    p_data_table         regclass, -- data table
    p_data_id_col        text,     -- identifier column
    p_data_feat_col      text,     -- feature vector array column
    p_omni_table         regclass, -- Omni-coordinate table
    p_omni_id_col        text,     -- Omni identifier column
    p_omni_array_col     text,     -- Omni-coordinate array column
    p_pivots_table       regclass, -- pivot table
    p_center_id          bigint,   -- query center id
    p_r_cap              float8,   -- initial query radius
    p_k                  int,      -- number of diversified neighbors
    p_m                  int       -- rings are indexed from 0 to m
)
returns table (
    id              bigint, -- returned object identifier
    dist_center     float8, -- distance from the returned object to the query center
    pivot_radical   int[]   -- pivot-ring radical of the returned object
)
language plpgsql
as $$
declare
    -- Query center feature vector
    v_center_feat   float8[];

    -- Query center Omni-coordinates
    v_center_omni   float8[];

    -- Ring configuration
    v_max_dist      float8;
    v_band_size     float8;

    -- Dynamic query radius
    v_rq            float8;
    v_expand        float8;

    -- Number of accepted distinct radicals
    v_cnt           int := 0;

    -- Dynamic SQL used during expansion
    v_sql           text;
    rec             record;

    -- Auxiliary variables
    v_existing_dist float8;
    v_worst_rad     int[];
    v_rows_round    int;
begin
    if p_k is null or p_k <= 0 then
        raise exception 'p_k must be greater than 0.';
    end if;

    if p_m is null or p_m < 0 then
        raise exception 'p_m must be greater than or equal to 0.';
    end if;

    -- Load the query center feature vector
    execute format(
        'select %1$I from %2$s where %3$I = $1',
        p_data_feat_col,
        p_data_table,
        p_data_id_col
    )
    into v_center_feat
    using p_center_id;

    if v_center_feat is null then
        raise exception 'Center id % was not found in table %.', p_center_id, p_data_table;
    end if;

    -- Load the query center Omni-coordinates
    execute format(
        'select %1$I from %2$s where %3$I = $1',
        p_omni_array_col,
        p_omni_table,
        p_omni_id_col
    )
    into v_center_omni
    using p_center_id;

    if v_center_omni is null then
        raise exception 'Center id % was not found in Omni table %.', p_center_id, p_omni_table;
    end if;

    -- Load max_dist and compute the ring thickness
    select dpm.max_dist
      into v_max_dist
      from dataset_pivot_metadata dpm
     where dpm.data_table = p_data_table
       and dpm.pivots_table = p_pivots_table;

    if v_max_dist is null then
        raise exception
            'dataset_pivot_metadata does not contain max_dist for data_table=% and pivots_table=%.',
            p_data_table, p_pivots_table;
    end if;

    v_band_size := v_max_dist / (p_m + 1);

    if v_band_size <= 0 then
        raise exception
            'Invalid band_size: max_dist=% and (m+1)=%.',
            v_max_dist, (p_m + 1);
    end if;

     -- Temporary table storing the best representative per radical
    create temporary table if not exists tmp_best_rad (
        pivot_radical int[] primary key,
        id              bigint,
        dist_center     float8
    ) on commit drop;

    truncate tmp_best_rad;

    -- Temporary table storing visited object identifiers
    create temporary table if not exists tmp_seen_ids (
        id bigint primary key
    ) on commit drop;

    truncate tmp_seen_ids;

    -- Initialize the query radius
    if p_r_cap is null then
        v_rq := v_band_size;
    else
        v_rq := p_r_cap;
    end if;

    v_expand := v_band_size;

    -- Expand the radius until k distinct radicals are found
    loop
        v_rows_round := 0;

        v_sql := format($SQL$
with
params as (
    select
        %1$L::float8[] as center_feat,
        %2$L::float8[] as center_omni,
        %3$L::float8   as r,
        %4$L::float8   as band_size,
        %5$L::int      as m
),

-- Omni pruning
candidates as (
    select
        o.%7$I::bigint as id,
        o.%8$I         as omni_coords
    from %6$s o
    cross join params p
    where o.%7$I <> %9$L
      and not exists (
            select 1
            from tmp_seen_ids s
            where s.id = o.%7$I
      )
      and (
            select bool_and(
                abs(o.%8$I[i] - p.center_omni[i]) <= p.r
            )
            from generate_subscripts(o.%8$I, 1) g(i)
          )
),

-- Real distance computation
real_dist as (
    select
        d.%10$I::bigint as id,
        s.dist_center,
        c.omni_coords
    from %11$s d
    join candidates c
      on c.id = d.%10$I
    cross join params p
    join lateral (
        select sqrt(sum((a - b) * (a - b)))::float8 as dist_center
        from unnest(d.%12$I, p.center_feat) as t(a, b)
    ) s on true
    where s.dist_center <= p.r
),

-- Build pivot-ring radicals
radicals as (
    select
        id,
        dist_center,
        array(
            select least(
                floor(real_dist.omni_coords[i] / p.band_size)::int,
                p.m
            )
            from generate_subscripts(real_dist.omni_coords, 1) g(i)
            cross join params p
            order by i
        ) as pivot_radical
    from real_dist
)

select
    id,
    dist_center,
    pivot_radical
from radicals
order by dist_center, id
$SQL$,
            v_center_feat,     -- %1$L
            v_center_omni,     -- %2$L
            v_rq,              -- %3$L
            v_band_size,       -- %4$L
            p_m,               -- %5$L
            p_omni_table,      -- %6$s
            p_omni_id_col,     -- %7$I
            p_omni_array_col,  -- %8$I
            p_center_id,       -- %9$L
            p_data_id_col,     -- %10$I
            p_data_table,      -- %11$s
            p_data_feat_col    -- %12$I
        );

        for rec in execute v_sql loop
            v_rows_round := v_rows_round + 1;

            -- Mark the object as visited
            insert into tmp_seen_ids(id)
            values (rec.id)
            on conflict do nothing;

             -- Keep the closest object for each radical
            select b.dist_center
              into v_existing_dist
              from tmp_best_rad b
             where b.pivot_radical = rec.pivot_radical;

            if v_existing_dist is null then
                insert into tmp_best_rad(pivot_radical, id, dist_center)
                values (rec.pivot_radical, rec.id, rec.dist_center)
                on conflict do nothing;

                if found then
                    v_cnt := v_cnt + 1;

                    -- Remove the farthest radical if more than k are stored
                    if v_cnt > p_k then
                        select b.pivot_radical
                          into v_worst_rad
                          from tmp_best_rad b
                         order by b.dist_center desc, b.id desc
                         limit 1;

                        delete from tmp_best_rad b
                         where b.pivot_radical = v_worst_rad;

                        v_cnt := p_k;
                    end if;
                end if;
            end if;
        end loop;

        --- Stop when k radicals have been collected
        if v_cnt = p_k then
            select max(b.dist_center)
              into v_rq
              from tmp_best_rad b;
            exit;
        end if;

        -- Expand the query radius
        v_rq := v_rq + v_expand;

        if v_rq > v_max_dist then
            v_rq := v_max_dist;
        end if;

        -- Stop if no new objects were found at the maximum radius
        if v_rows_round = 0 and v_rq >= v_max_dist then
            exit;
        end if;

        if v_rq >= v_max_dist and v_rows_round = 0 then
            exit;
        end if;
    end loop;

    -- Return the selected representatives
    return query
    select
        b.id,
        b.dist_center,
        b.pivot_radical
    from tmp_best_rad b
    order by b.dist_center, b.id
    limit p_k;

end;
$$;
