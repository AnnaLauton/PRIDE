create or replace function pride_r_seq(
    p_table              regclass, -- data table
    p_id_col             text,     -- identifier column
    p_feat_array_col     text,     -- feature vector array column
    p_pivots_table       regclass, -- pivot table
    p_pivot_pos_col      text,     -- pivot ordering column
    p_pivot_array_col    text,     -- pivot feature vector array column
    p_center_id          bigint,   -- query center id
    p_r                  float8,   -- query radius
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
    -- Query center
    v_center          float8[];

    -- Ring configuration
    v_max_dist        float8;
    v_band_size       float8;

    -- Number of pivots
    v_num_pivots      int;

    -- Dynamic SQL for the candidate scan
    v_sql_scan        text;
    rec               record;

    -- Radical construction
    v_rad             int[];
    v_ring            int;
    v_dist_cp         float8;
    v_pivot_feat      float8[];
    i                 int;
begin
    if p_r is null or p_r < 0 then
        raise exception 'p_r must be greater than or equal to 0.';
    end if;

    if p_m is null or p_m < 0 then
        raise exception 'p_m must be greater than or equal to 0.';
    end if;

    -- Load the query center
    execute format(
        'select %1$I from %2$s where %3$I = $1',
        p_feat_array_col,
        p_table,
        p_id_col
    )
    into v_center
    using p_center_id;

    if v_center is null then
        raise exception 'center_id % was not found in table %.', p_center_id, p_table;
    end if;

    -- Load max_dist and compute the ring thickness
    select dpm.max_dist
      into v_max_dist
      from dataset_pivot_metadata dpm
     where dpm.data_table = p_table
       and dpm.pivots_table = p_pivots_table;

    if v_max_dist is null then
        raise exception
            'dataset_pivot_metadata does not contain max_dist for data_table=% and pivots_table=%.',
            p_table, p_pivots_table;
    end if;

    v_band_size := (v_max_dist / (p_m + 1))::float8;

    if v_band_size <= 0 then
        raise exception
            'Invalid band_size: max_dist=% and (m+1)=%.',
            v_max_dist, (p_m + 1);
    end if;

    -- Load pivots into a temporary table
    create temporary table if not exists tmp_pivot_arrays (
        pivot_pos   int primary key,
        feat_array  float8[] not null
    ) on commit drop;

    truncate tmp_pivot_arrays;

    execute format($fmt$
        insert into tmp_pivot_arrays(pivot_pos, feat_array)
        select
            p.%1$I::int,
            p.%2$I::float8[]
        from %3$s p
        order by p.%1$I
    $fmt$,
        p_pivot_pos_col,
        p_pivot_array_col,
        p_pivots_table
    );

    select count(*)
      into v_num_pivots
      from tmp_pivot_arrays;

    if v_num_pivots < 1 then
        raise exception 'No pivots were found in table %.', p_pivots_table;
    end if;

    -- Temporary table storing the best representative per radical
    create temporary table if not exists tmp_best_rad (
        pivot_radical	int[] primary key,
        id              bigint,
        dist_center     float8
    ) on commit drop;

    truncate tmp_best_rad;

    -- Scan candidates within the query radius
    v_sql_scan := format($fmt$
        with q as (
            select
                %1$L::bigint  as center_id,
                %5$L::float8  as r,
                $1::float8[]  as center_feat
        )
        select
            d.%2$I::bigint as id,
            d.%3$I         as feat,
            dist_l2(q.center_feat, d.%3$I)::float8 as dist_center
        from %4$s d
        cross join q
        where d.%2$I <> q.center_id
          and dist_l2(q.center_feat, d.%3$I) <= q.r
        order by dist_l2(q.center_feat, d.%3$I), d.%2$I
    $fmt$,
        p_center_id,
        p_id_col,
        p_feat_array_col,
        p_table,
        p_r
    );

    for rec in execute v_sql_scan using v_center loop
        -- Build the pivot-based radical
        v_rad := array[]::int[];

        for i in 1..v_num_pivots loop
            select t.feat_array
              into v_pivot_feat
              from tmp_pivot_arrays t
             where t.pivot_pos = i;

            v_dist_cp := dist_l2(rec.feat, v_pivot_feat);
            v_ring := least(floor(v_dist_cp / v_band_size)::int, p_m);
            v_rad := v_rad || v_ring;
        end loop;

        -- Keep the closest object for each radical
        insert into tmp_best_rad(pivot_radical, id, dist_center)
        values (v_rad, rec.id, rec.dist_center)
        on conflict do nothing;
    end loop;

    -- Return one representative per radical
    return query
    select
        b.id,
        b.dist_center,
        b.pivot_radical
    from tmp_best_rad b
    order by b.dist_center, b.id;

end;
$$;
