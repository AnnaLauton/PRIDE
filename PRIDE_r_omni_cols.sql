create or replace function pride_r_omni_cols(
    p_data_table         regclass, -- data table
    p_data_id_col        text,     -- identifier column
    p_data_feat_col      text,     -- feature vector array column
    p_omni_cols_table    regclass, -- column-based Omni-coordinate table
    p_omni_id_col        text,     -- Omni identifier column
    p_pivots_table       regclass, -- pivot table
    p_num_pivots         int,      -- number of pivot columns
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
	-- Query center feature vector
    v_center_feat      float8[];

	-- Query center Omni-coordinates
    v_center_omni      float8[];

	-- Ring configuration
    v_max_dist         float8;
    v_band_size        float8;

	-- Dynamic SQL used to load Omni-coordinates
    v_center_sql       text;

	-- Dynamic SQL used for query execution
    v_sql              text;

	-- Dynamic Omni predicates and radical expression
    v_where_pred       text := '';
    v_rad_expr         text := '';
    i                  int;
begin
    if p_r is null or p_r < 0 then
        raise exception 'p_r must be greater than or equal to 0.';
    end if;

    if p_m is null or p_m < 0 then
        raise exception 'p_m must be greater than or equal to 0.';
    end if;

    if p_num_pivots is null or p_num_pivots <= 0 then
        raise exception 'p_num_pivots must be greater than 0.';
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
    v_center_sql := 'select array[';

    for i in 1..p_num_pivots loop
        if i > 1 then
            v_center_sql := v_center_sql || ', ';
        end if;

        v_center_sql := v_center_sql || format('pivot_%s', i);
    end loop;

    v_center_sql := v_center_sql || format(
        ']::float8[] from %s where %I = $1',
        p_omni_cols_table,
        p_omni_id_col
    );

    execute v_center_sql
    into v_center_omni
    using p_center_id;

    if v_center_omni is null then
        raise exception 'Center id % was not found in Omni columns table %.', p_center_id, p_omni_cols_table;
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

    -- Build Omni pruning predicates and radical expression
    for i in 1..p_num_pivots loop
        if i > 1 then
            v_where_pred := v_where_pred || ' and ';
            v_rad_expr := v_rad_expr || ', ';
        end if;

        v_where_pred := v_where_pred || format(
            'o.pivot_%s between (%L::float8 - %L::float8) and (%L::float8 + %L::float8)',
            i, v_center_omni[i], p_r, v_center_omni[i], p_r
        );

        v_rad_expr := v_rad_expr || format(
            'least(floor(o.pivot_%s / %L::float8)::int, %s)',
            i, v_band_size, p_m
        );
    end loop;

    -- Main query
    v_sql := format($SQL$
with
params as (
    select
        %1$L::float8[] as center_feat,
        %2$L::float8   as r
),

-- Omni pruning using one predicate per pivot column
candidates as (
    select
        o.%4$I::bigint as id,
        array[%5$s]::int[] as pivot_radical
    from %3$s o
    where o.%4$I <> %6$L
      and %7$s
),

-- Real distance computation
real_dist as (
    select
        d.%8$I::bigint as id,
        s.dist_center,
        c.pivot_radical
    from %9$s d
    join candidates c
      on c.id = d.%8$I
    cross join params p
    join lateral (
        select sqrt(sum((a - b) * (a - b)))::float8 as dist_center
        from unnest(d.%10$I, p.center_feat) as t(a, b)
    ) s on true
    where s.dist_center <= p.r
),

-- Keep the closest object for each radical
ranked as (
    select
        id,
        dist_center,
        pivot_radical,
        row_number() over (
            partition by pivot_radical
            order by dist_center, id
        ) as rn
    from real_dist
)

select
    id::bigint,
    dist_center,
    pivot_radical
from ranked
where rn = 1
order by dist_center, id
$SQL$,
        v_center_feat,    -- %1$L
        p_r,              -- %2$L
        p_omni_cols_table,-- %3$s
        p_omni_id_col,    -- %4$I
        v_rad_expr,       -- %5$s
        p_center_id,      -- %6$L
        v_where_pred,     -- %7$s
        p_data_id_col,    -- %8$I
        p_data_table,     -- %9$s
        p_data_feat_col   -- %10$I
    );

    return query execute v_sql;
end;
$$;
