create or replace function pride_k_omni_cols(
    p_data_table         regclass, -- data table
    p_data_id_col        text,     -- identifier column
    p_data_feat_col      text,     -- feature vector array column
    p_omni_cols_table    regclass, -- column-based Omni-coordinate table
    p_omni_id_col        text,     -- Omni identifier column
    p_pivots_table       regclass, -- pivot table
    p_num_pivots         int,      -- number of pivot columns
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
    v_center_feat      float8[];

     -- Query center Omni-coordinates
    v_center_omni      float8[];

    -- Ring configuration
    v_max_dist         float8;
    v_band_size        float8;

    -- Dynamic query radius
    v_rq               float8;
    v_expand           float8;

    -- Number of accepted distinct radicals
    v_cnt              int := 0;

    -- Dynamic SQL fragments
    v_center_sql       text;
    v_where_pred       text;
    v_rad_expr         text;
    v_sql              text;

    -- Loop variables
    i                  int;
    rec                record;

    -- Auxiliary variables
    v_existing_dist    float8;
    v_worst_rad        int[];
    v_rows_round       int;
begin
    if p_k is null or p_k <= 0 then
        raise exception 'p_k must be greater than 0.';
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

    -- Initialize the query radius
    if p_r_cap is null then
        v_rq := v_band_size;
    else
        v_rq := p_r_cap;
    end if;

    v_expand := v_band_size;

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

    -- Expand the radius until k distinct radicals are found
    loop
        v_rows_round := 0;
        v_where_pred := '';
        v_rad_expr := '';

        -- Build Omni pruning predicates and radical expression
        for i in 1..p_num_pivots loop
            if i > 1 then
                v_where_pred := v_where_pred || ' and ';
                v_rad_expr := v_rad_expr || ', ';
            end if;

            v_where_pred := v_where_pred || format(
                'o.pivot_%s between (%L::float8 - %L::float8) and (%L::float8 + %L::float8)',
                i, v_center_omni[i], v_rq, v_center_omni[i], v_rq
            );

            v_rad_expr := v_rad_expr || format(
                'least(floor(o.pivot_%s / %L::float8)::int, %s)',
                i, v_band_size, p_m
            );
        end loop;

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
      and not exists (
            select 1
            from tmp_seen_ids s
            where s.id = o.%4$I
      )
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
)

select
    id,
    dist_center,
    pivot_radical
from real_dist
order by dist_center, id
$SQL$,
            v_center_feat,     -- %1$L
            v_rq,              -- %2$L
            p_omni_cols_table, -- %3$s
            p_omni_id_col,     -- %4$I
            v_rad_expr,        -- %5$s
            p_center_id,       -- %6$L
            v_where_pred,      -- %7$s
            p_data_id_col,     -- %8$I
            p_data_table,      -- %9$s
            p_data_feat_col    -- %10$I
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

        -- Stop when k radicals have been collected
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
